#!/usr/bin/env python
import cdsapi
c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete', {
    'class': 'ea',
    'date': '2021-12-01/to/2022-01-31',
    'expver': '1',
    'levelist': '2000',
    'levtype': 'pv',
    'param': '3.128/54.128/131.128/132.128',
    'stream': 'oper',
    'time': '00:00:00/06:00:00/12:00:00/18:00:00',
    'type': 'an',
    'grid': [0.25,0.25],
    'format': 'netcdf',
}, 'era5_2pvu_2021120100_2022013118.nc')


