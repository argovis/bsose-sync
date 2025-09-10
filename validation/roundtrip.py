from pymongo import MongoClient
import xarray, numpy

client = MongoClient('mongodb://database/argo')
db = client.argo

# open some file objects for the upstream netcdf files
theta13 = xarray.open_dataset('data/Theta_bsoseI156_2013_5day.nc')
theta14 = xarray.open_dataset('data/Theta_bsoseI156_2014_5day.nc')
salt13 = xarray.open_dataset('data/Salt_bsoseI156_2013_5day.nc')
salt14 = xarray.open_dataset('data/Salt_bsoseI156_2014_5day.nc')

longitudes = [float(x) for x in list(theta13.coords['XC'].values)]
latitudes = [float(x) for x in list(theta13.coords['YC'].values)]
levels = [float(x) for x in list(theta13.coords['Z'].values)]

while True:
    logmessage = ''
    
    # get a random profile, or pick one by ID
    p = list(db.bsose.aggregate([{"$sample": {"size": 1}}]))[0]
    #p = list(db.bsose.find({"_id":"xxx"}))[0]
    m = list(db.timeseriesMeta.find({"_id":p['metadata'][0]}))[0]
    logmessage += 'Checking profile id ' + str(p['_id']) + '\n'
    loglen = len(logmessage)

    # extract some metadata
    lat = p['geolocation']['coordinates'][1]
    lon = p['geolocation']['coordinates'][0]%360
    level = p['level']

    # find variable indexes
    try:
        theta_idx = p['data_info'][0].index('THETA')
    except ValueError:
        theta_idx = -1
    try:
        salt_idx = p['data_info'][0].index('SALT')
    except ValueError:
        salt_idx = -1
    
    # extract variable vectors
    if theta_idx > -1:
        theta = p['data'][theta_idx]
    else:
        theta = []
    if salt_idx > -1:
        salt = p['data'][salt_idx]
    else:
        salt = []

    # find this subset in the netcdf files
    lonidx = longitudes.index(lon)
    latidx = latitudes.index(lat)
    levelidx = levels.index(-1*level)

    # extract timeseries at these indexes
    theta_timeseries = list(theta13['THETA'][:,levelidx,latidx,lonidx].values) + list(theta14['THETA'][:,levelidx,latidx,lonidx].values)
    theta_timeseries = [float(x) for x in theta_timeseries]
    salt_timeseries = list(salt13['SALT'][:,levelidx,latidx,lonidx].values) + list(salt14['SALT'][:,levelidx,latidx,lonidx].values)
    salt_timeseries = [float(x) for x in salt_timeseries]

    if not numpy.allclose(theta, theta_timeseries, equal_nan=True):
        logmessage += 'Theta mismatch\n'
        logmessage += 'Profile ID: ' + str(p['_id']) + '\n'
    if not numpy.allclose(salt, salt_timeseries, equal_nan=True):
        logmessage += 'Salt mismatch\n'
        logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

    for f in [theta13, theta14, salt13, salt14]:
         # check metadata doc against all upstream files - assumed some things are a function only of lon,lat, recorded only once

        if f['rA'].values[latidx,lonidx] != m['cell_area']:
            logmessage += 'Cell area mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'
    
        if f['Depth'].values[latidx,lonidx] != m['ocean_depth']:
            logmessage += 'ocean_depth mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        if f['maskInC'].values[latidx,lonidx] != m['interior_2d_mask']:
            logmessage += 'interior_2d_mask mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        if f['rSurfC'].values[latidx, lonidx] != m['depth_r0_to_ref_surface']:
            logmessage += 'depth_r0_to_ref_surface mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        if f['rLowC'].values[latidx, lonidx] != m['depth_r0_to_bottom']:
            logmessage += 'depth_r0_to_bottom mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        # same idea for data doc
        if f['drF'].values[levelidx] != p['cell_z_size']:
            logmessage += 'cell_z_size mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        if f['hFacC'].values[levelidx,latidx,lonidx] != p['cell_vertical_fraction']:
            logmessage += 'cell_vertical_fraction mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'
        
        if f['maskC'].values[levelidx,latidx,lonidx] != p['sea_binary_mask_at_t_locaiton']:
            logmessage += 'sea_binary_mask_at_t_locaiton mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

        if f['rhoRef'].values[levelidx] != p['reference_density_profile']:
            logmessage += 'reference_density_profile mismatch\n'
            logmessage += 'Profile ID: ' + str(p['_id']) + '\n'

    if len(logmessage) > loglen:
        print(logmessage)

    
