import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
            'mean_sea_level_pressure', 'surface_net_solar_radiation', 'surface_net_solar_radiation_clear_sky',
            'surface_net_thermal_radiation', 'surface_net_thermal_radiation_clear_sky', 'surface_solar_radiation_downward_clear_sky',
            'surface_solar_radiation_downwards', 'surface_thermal_radiation_downward_clear_sky', 'surface_thermal_radiation_downwards',
            'total_cloud_cover', 'total_column_water_vapour', 'vertical_integral_of_northward_water_vapour_flux',
        ],
        'year': '1983',
        'month': [
        '04', '05',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '06:00', '12:00',
            '18:00',
        ],
        'format': 'netcdf',
    },
    'era5_sfclevel_1983040100_1983053018.nc')