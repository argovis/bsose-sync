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
use mongodb::bson::Bson;

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

    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let dv = &args[2];
    let lolat = args[3].parse::<usize>()?;
    let hilat = args[4].parse::<usize>()?;
    let lolong = args[5].parse::<usize>()?;
    let hilong = args[6].parse::<usize>()?;

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

    // Rust structs to describe documents in the "bsose" collections
    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct Sourcedoc {
        source: Vec<String>,
        file: String
    }

    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct BsoseMetadoc {
        _id: String,
        data_type: String,
        date_updated_argovis: DateTime,
        timeseries: Vec<DateTime>,
        source: Vec<Sourcedoc>,
        cell_area: f64,
        ocean_depth: f64,
        depth_r0_to_bottom: f64,
        interior_2d_mask: bool,
        depth_r0_to_ref_surface: f64
    }

    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct BsoseDocument {
        _id: String,
        metadata: Vec<String>,
        basin: i32,
        geolocation: Geolocation,
        level: f64,
        data: Vec<Vec<f64>>,
        data_info: (Vec<String>, Vec<String>, Vec<Vec<String>>),
        cell_vertical_fraction: f64,
        sea_binary_mask_at_t_locaiton: bool,
        ctrl_vector_3d_mask: bool,
        cell_z_size: f64,
        reference_density_profile: f64,
    }

    #[derive(Serialize, Deserialize, Debug, Clone)]
    struct Geolocation {
        #[serde(rename = "type")]
        location_type: String,
        coordinates: [f64; 2],
    }

    // collection objects
    let bsose = client.database("argo").collection::<BsoseDocument>("bsoseX");
    //let bsose_meta = client.database("argo").collection::<BsoseMetadoc>("timeseriesMeta");
  
    let file = netcdf::open(filename)?;

    // basin lookup
    let basinfile = netcdf::open("/tmp/basinmask_01.nc")?;
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
    let datavar = &file.variable(dv).expect("Could not find data variable");
    let mut units: String = String::from("");
    let mut long_name: String = String::from("");
    if let netcdf::AttrValue::Str(u) = datavar.attribute_value("units").unwrap()? {
        units = u;
    }
    if let netcdf::AttrValue::Str(u) = datavar.attribute_value("long_name").unwrap()? {
        long_name = u;
    }

    // construct metadata
    let n_timesteps = time.len();
    let mut timeseries = Vec::new();
    for timeidx in 0..n_timesteps {
        timeseries.push(bson::DateTime::parse_rfc3339_str((t0 + Duration::seconds(time.value::<i64, _>(timeidx)?)).to_rfc3339().replace("+00:00", "Z")).unwrap());
    }
/*
    for latidx in lolat..hilat {
        for lonidx in lolong..hilong {
            let lon_val = tidylon(lon.value::<f64, _>([lonidx])?);

            // construct metadata documents
            let metaid = format!("{:.3}_{:.3}", lon_val, lat.value::<f64, _>([latidx])?);
            bsose_meta.insert_one(BsoseMetadoc{
                _id: metaid.clone(),
                data_type: String::from("BSOSE-profile"),
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
            }, None).await?;

        }
    }
*/

    for latidx in lolat..hilat { //588 //lat.len() {
        let lat_val = lat.value::<f64, _>([latidx])?;
        for lonidx in lolong..hilong {
            let lon_val = tidylon(lon.value::<f64, _>([lonidx])?);
            // construct data documents, one timeseries per lon/lat/level triple
            for levelidx in 0..depth.len() {
                let basin = find_basin(&basins, lon_val, lat_val);
                let mut datavar_profile = Vec::new();
                for timeidx in 0..n_timesteps {
                    datavar_profile.push(datavar.value::<f64, _>([timeidx, levelidx, latidx, lonidx])? as f64);
                }
                let id = format!("{:.3}_{:.3}_{:.3}", lon_val, lat_val, depth.value::<f64, _>(levelidx)?);

                // Check if a document with property "_id" matching id exists
                let existing_doc = bsose.find_one(doc! { "_id": id.clone() }, None).await?;

                if let Some(mut doc) = existing_doc {
                    // Append the value of datavar_profile to the existing "data" property
                    doc.data.push(datavar_profile.clone());
                    doc.data_info.0.push(dv.to_string());
                    doc.data_info.2.push(vec!(units.clone(), long_name.clone()));
                    let filter = doc! {"_id": id };
                    bsose.replace_one(filter, doc, None).await?;
                } else {
                    if !datavar_profile.iter().all(|&x| x == 0.0) {
                        bsose.insert_one(BsoseDocument {
                            _id: id,
                            metadata: vec![format!("{:.3}_{:.3}", lon_val, lat_val)],
                            basin: basin,
                            geolocation: Geolocation{
                                location_type: String::from("Point"),
                                coordinates: [lon_val, lat_val]
                            },
                            level: -1.0 * depth.value::<f64, _>(levelidx)?,
                            data: vec![datavar_profile.clone()],
                            data_info: (
                                vec!(dv.to_string()), 
                                vec!(String::from("units"), String::from("long_name")),
                                vec!(
                                    vec!(units.clone(), long_name.clone())
                                )
                            ),
                            cell_vertical_fraction: cell_vertical_fraction.value::<f64, _>((levelidx, latidx, lonidx))?,
                            sea_binary_mask_at_t_locaiton: sea_binary_mask_at_t_locaiton.value::<i8, _>((levelidx, latidx, lonidx))? != 0,
                            ctrl_vector_3d_mask:  ctrl_vector_3d_mask.value::<i8, _>((levelidx, latidx, lonidx))? != 0,
                            cell_z_size: cell_z_size.value::<f64, _>(levelidx)?,
                            reference_density_profile: reference_density_profile.value::<f64, _>(levelidx)?
                        }, None).await?;
                    }
                }
            }
        }
    }

    Ok(())
}
