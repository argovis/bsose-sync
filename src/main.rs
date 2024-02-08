use netcdf;
use tokio;
use std::error::Error;
use chrono::Utc;
use chrono::TimeZone;
use chrono::Duration;
use std::env;
use mongodb::bson::{doc};
use mongodb::bson::DateTime;
use mongodb::{Client, options::{ClientOptions, ResolverConfig}};
use serde::{Deserialize, Serialize};

fn tidylon(longitude: f64) -> f64{
    // map longitude on [0,360] to [-180,180], required for mongo indexing
    if longitude <= 180.0{
        return longitude;
    }
    else{
        return longitude-360.0;
    }
}

fn find_basin(basins: &netcdf::Variable, longitude: f64, latitude: f64) -> i32 {    
    let lonplus = (longitude-0.5).ceil()+0.5;
    let lonminus = (longitude-0.5).floor()+0.5;
    let latplus = (latitude-0.5).ceil()+0.5;
    let latminus = (latitude-0.5).floor()+0.5;

    let lonplus_idx = (lonplus - -179.5) as usize;
    let lonminus_idx = (lonminus - -179.5) as usize;
    let latplus_idx = (latplus - -77.5) as usize;
    let latminus_idx = (latminus - -77.5) as usize;

    let corners_idx = [
        // bottom left corner, clockwise
        [latminus_idx, lonminus_idx],
        [latplus_idx, lonminus_idx],
        [latplus_idx, lonplus_idx],
        [latminus_idx, lonplus_idx]
    ];

    let distances = [
        (f64::powi(longitude-lonminus, 2) + f64::powi(latitude-latminus, 2)).sqrt(),
        (f64::powi(longitude-lonminus, 2) + f64::powi(latitude-latplus, 2)).sqrt(),
        (f64::powi(longitude-lonplus, 2) + f64::powi(latitude-latplus, 2)).sqrt(),
        (f64::powi(longitude-lonplus, 2) + f64::powi(latitude-latminus, 2)).sqrt()
    ];

    let mut closecorner_idx = corners_idx[0];
    let mut closedist = distances[0];
    for i in 1..4 {
        if distances[i] < closedist{
            closecorner_idx = corners_idx[i];
            closedist = distances[i];
        }
    }

    match basins.value::<i64,_>(closecorner_idx){
        Ok(idx) => idx as i32,
        Err(e) => panic!("basin problems: {:?} {:#?}", e, closecorner_idx)
    }   
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {

    // setup /////////////////////////////////////////////////

    // mongodb setup
    // Load the MongoDB connection string from an environment variable:
    let client_uri =
       env::var("MONGODB_URI").expect("You must set the MONGODB_URI environment var!"); 

    // A Client is needed to connect to MongoDB:
    // An extra line of code to work around a DNS issue on Windows:
    let options =
       ClientOptions::parse_with_resolver_config(&client_uri, ResolverConfig::cloudflare())
          .await?;
    let client = Client::with_options(options)?; 

    // collection objects
    let bsose = client.database("argo").collection("bsoseX");
    let bsose_meta = client.database("argo").collection("bsoseMetaX");
  
    // Rust structs to serialize time properly
    #[derive(Serialize, Deserialize, Debug)]
    struct Sourcedoc {
        source: Vec<String>,
        file: String
    }

    #[derive(Serialize, Deserialize, Debug)]
    struct BsoseMetadoc {
        _id: String,
        data_type: String,
        data_info: (Vec<String>, Vec<String>, Vec<Vec<String>>),
        date_updated_argovis: DateTime,
        timeseries: Vec<DateTime>,
        source: Vec<Sourcedoc>,
        cell_area: f64,
        ocean_depth: f64,
        depth_r0_to_bottom: f64,
        interior_2d_mask: bool,
        depth_r0_to_ref_surface: f64
    }


    let file = netcdf::open("data/O2_bsoseI139_2013to2021_5dy.nc")?;

    // basin lookup
    let basinfile = netcdf::open("data/basinmask_01.nc")?;
    let basins = &basinfile.variable("BASIN_TAG").expect("Could not find variable 'BASIN_TAG'");

    // all times recorded as days since Dec 1 2012
    let t0 = Utc.with_ymd_and_hms(2012, 12, 1, 0, 0, 0).unwrap();

    // document construction //////////////////////////////////////

    // variable extraction
    let lat = &file.variable("YC").expect("Could not find variable 'YC'");
    let lon = &file.variable("XC").expect("Could not find variable 'XC'");
    let depth = &file.variable("Z").expect("Could not find variable 'Z'");
    let time = &file.variable("time").expect("Could not find variable 'time'");
    let cell_area = &file.variable("rA").expect("Could not find variable 'rA'");
    let ocean_depth = &file.variable("Depth").expect("Could not find variable 'Depth'");
    let depth_r0_to_bottom = &file.variable("rLowC").expect("Could not find variable 'rLowC'");
    let interior_2d_mask = &file.variable("maskInC").expect("Could not find variable 'maskInC'");
    let depth_r0_to_ref_surface = &file.variable("rSurfC").expect("Could not find variable 'rSurfC'");
    let cell_vertical_fraction = &file.variable("hFacC").expect("Could not find variable 'hFacC'");
    let sea_binary_mask_at_t_locaiton = &file.variable("maskC").expect("Could not find variable 'maskC'");
    let ctrl_vector_3d_mask = &file.variable("maskCtrlC").expect("Could not find variable 'maskCtrlC'");
    let cell_z_size = &file.variable("drF").expect("Could not find variable 'drF'");
    let reference_density_profile = &file.variable("rhoRef").expect("Could not find variable 'rhoRef'");
    let track03 = &file.variable("TRAC03").expect("Could not find variable 'TRAC03'");

    // construct metadata
    let n_timesteps = time.len();
    let mut timeseries = Vec::new();
    for timeidx in 0..n_timesteps {
        timeseries.push(bson::DateTime::parse_rfc3339_str((t0 + Duration::seconds(time.value::<i64, _>(timeidx)?)).to_rfc3339().replace("+00:00", "Z")).unwrap());
    }

    for latidx in 0..lat.len() {
        for lonidx in 0..lon.len() {
            let lon_val = tidylon(lon.value::<f64, _>([lonidx])?);

            // construct metadata documents
            let metaid = format!("{:.3}_{:.3}", lon_val, lat.value::<f64, _>(latidx)?);

            let metadata = BsoseMetadoc {
                _id: metaid.clone(),
                data_type: String::from("BSOSE-profile"),
                data_info: (
                    vec!(String::from("TRAC03")), 
                    vec!(String::from("units"), String::from("long_name")),
                    vec!(
                        vec!(String::from("mol O/m"), String::from("Dissolved Oxygen concentration"))
                    )
                ),
                date_updated_argovis: DateTime::now(),
                timeseries: timeseries.clone(),
                source: vec!(
                    Sourcedoc{
                        source: vec!(String::from("BSOSE")),
                        file: String::from("O2_bsoseI139_2013to2021_5dy.nc")
                    }
                ),
                cell_area: cell_area.value::<f64, _>((latidx, lonidx))?,
                ocean_depth: ocean_depth.value::<f64, _>((latidx, lonidx))?,
                depth_r0_to_bottom: depth_r0_to_bottom.value::<f64, _>((latidx, lonidx))?,
                interior_2d_mask: interior_2d_mask.value::<i8, _>((latidx, lonidx))? != 0,
                depth_r0_to_ref_surface: depth_r0_to_ref_surface.value::<f64, _>((latidx, lonidx))?,
            };

            //println!("{:?}", metadata);
            let metadata_doc = bson::to_document(&metadata).unwrap();
            bsose_meta.insert_one(metadata_doc.clone(), None).await?;
        }
    }

    for latidx in 0..lat.len() {
        let lat_val = lat.value::<f64, _>([latidx])?;
        let mut docs = Vec::new(); // collect all the docs for this latitude, and write all at once.
        for lonidx in 0..lon.len() {
            let lon_val = tidylon(lon.value::<f64, _>([lonidx])?);
            // construct data documents, one timeseries per lon/lat/level triple
            for levelidx in 0..depth.len() {
                let basin = find_basin(&basins, lon_val, lat_val);
                let mut track03_profile = Vec::new();
                for timeidx in 0..n_timesteps {
                    track03_profile.push(track03.value::<f64, _>([timeidx, levelidx, latidx, lonidx])? as f64);
                }
                let data = doc!{
                    "_id": format!("{:.3}_{:.3}_{:.3}", lon_val, lat_val, depth.value::<f64, _>(levelidx)?),
                    "metadata": [format!("{:.3}_{:.3}", lon_val, lat_val)],
                    "basin": basin,
                    "geolocation": {
                        "type": "Point",
                        "coordinates": [lon_val, lat_val]
                    },
                    "level": depth.value::<f64, _>(levelidx)?,
                    "data": [track03_profile.clone()],
                    "cell_vertical_fraction": cell_vertical_fraction.value::<f64, _>((levelidx, latidx, lonidx))?,
                    "sea_binary_mask_at_t_locaiton": sea_binary_mask_at_t_locaiton.value::<i8, _>((levelidx, latidx, lonidx))? != 0,
                    "ctrl_vector_3d_mask": ctrl_vector_3d_mask.value::<i8, _>((levelidx, latidx, lonidx))? != 0,
                    "cell_z_size": cell_z_size.value::<f64, _>(levelidx)?,
                    "reference_density_profile": reference_density_profile.value::<f64, _>(levelidx)?
                };
                docs.push(data);
            }
        }
        bsose.insert_many(docs, None).await?;
    }

    Ok(())
}
