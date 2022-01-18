#!/usr/bin/env python3

#
# Program to visualize the shockwave from Tonga using data provided by
# Japan Meteorological Agency.
#
# Copyright (C) 2022 Takashi Nakamoto
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

import numpy as np
from datetime import datetime
from math import nan, isnan
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfea
from pathlib import Path

class PressureData:
    def __init__(self, latitude_deg, longitude_deg, date_time, pressure_hPa, pressure_hPa_diff):
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.date_time     = date_time
        self.pressure_hPa  = pressure_hPa
        self.pressure_hPa_diff = pressure_hPa_diff

def read_pressure_data(filename):
    f = open( filename, 'r', encoding = 'utf-8' )
    
    lines = f.readlines()
    latitude_deg  = float( lines[4].split(',')[1] )
    longitude_deg = float( lines[5].split(',')[1] )

    data = []
    pressure_hPa_prev = nan

    for line in lines[9:]:
        fields = [r.strip() for r in line.split(',')]
        date_time    = datetime.strptime( fields[0], '%Y-%m-%d %H:%M' )
        pressure_hPa = float( fields[1] )
        pressure_hPa_diff = pressure_hPa - pressure_hPa_prev
        pressure_hPa_prev = pressure_hPa

        data.append( PressureData( latitude_deg, longitude_deg, date_time, pressure_hPa, pressure_hPa_diff ) )
    
    f.close()

    return data

# ==============================================================================
# Read presure data
# ==============================================================================
DATA_INPUT_DIR = Path( 'data' )
csv_files = DATA_INPUT_DIR.glob( '*.csv' )

records = []

for csv_file in csv_files:
    records.extend( read_pressure_data( csv_file ) )

# ==============================================================================
# Draw Japan
# ==============================================================================
fig = plt.figure( figsize = (16, 8) )
ax = fig.add_subplot( 1, 1, 1,
                      projection = ccrs.PlateCarree( central_longitude = 180 ) )
START_LATITUDE_DEG  =  24.0
END_LATITUDE_DEG    =  46.0
START_LONGITUDE_DEG = 120.0
END_LONGITUDE_DEG   = 160.0
ax.set_extent( ( START_LONGITUDE_DEG, END_LONGITUDE_DEG,
                 START_LATITUDE_DEG, END_LATITUDE_DEG ),
               ccrs.PlateCarree() )

ax.add_feature( cfea.OCEAN, color = '#00FFFF' )
ax.add_feature( cfea.LAND,  color = '#32CD32' )

# ==============================================================================
# Plot pressure data every 10 minutes
# ==============================================================================
FIG_OUTPUT_DIR = Path( 'figure' )
FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

MAX_PRESSURE_HPA_DIFF = 0.5

date_times = sorted( list( set( [ record.date_time for record in records ] ) ) )
for date_time in date_times:
    matched_records = [ record for record in records if record.date_time == date_time ]
    points = []
    for matched_record in matched_records:
        
        if isnan( matched_record.pressure_hPa_diff ):
            continue

        color = cm.rainbow( matched_record.pressure_hPa_diff / MAX_PRESSURE_HPA_DIFF / 2.0 + 0.5 )

        p = ax.plot( matched_record.longitude_deg,
                     matched_record.latitude_deg,
                     'o',
                     transform = ccrs.PlateCarree(),
                     markersize = 7,
                     color = color )
        points.append( p )

    date_time_str = date_time.strftime('%Y%m%d_%H%M')
    fig_output_filename = FIG_OUTPUT_DIR / f'{date_time_str}.png'
    fig.savefig( fig_output_filename )
    print( f'Generated {fig_output_filename}' )

    for p in points:
        pp = p.pop(0)
        pp.remove()
        del pp




