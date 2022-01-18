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
import matplotlib.lines as mlines
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
plt.rcParams['font.family'] = 'monospace'
fig = plt.figure( figsize = (10, 5) )
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
color_map = lambda dp: cm.rainbow( dp / MAX_PRESSURE_HPA_DIFF / 2.0 + 0.5 )

date_times = sorted( list( set( [ record.date_time for record in records ] ) ) )
for date_time in date_times[100:]:
    matched_records = [ record for record in records if record.date_time == date_time ]
    points = []
    longitude_degs = []
    latitude_degs  = []
    pressure_hPa_diffs = []
    
    for matched_record in matched_records:
        
        if isnan( matched_record.pressure_hPa_diff ):
            continue

        p = ax.plot( matched_record.longitude_deg,
                     matched_record.latitude_deg,
                     'o',
                     transform = ccrs.PlateCarree(),
                     markersize = 5,
                     color = color_map( matched_record.pressure_hPa_diff ) )
        points.append( p )

        if not isnan( matched_record.pressure_hPa_diff ):
            longitude_degs.append( matched_record.longitude_deg )
            latitude_degs.append( matched_record.latitude_deg )
            pressure_hPa_diffs.append( matched_record.pressure_hPa_diff )

# Removed contour plots because they are just ugly...
#
#    lines = ax.tricontour( longitude_degs, latitude_degs, pressure_hPa_diffs,
#                           linestyles = '-',
#                           colors = 'black',
#                           linewidth = 0.5,
#                           transform = ccrs.PlateCarree(),
#                           levels = np.linspace(-MAX_PRESSURE_HPA_DIFF * 1.1, MAX_PRESSURE_HPA_DIFF * 1.1, 23) )

    legend_items = []

    legend_title_lines = [f'{date_time.strftime("%Y-%m-%d %H:%M")} (JST)',
                           'Barometric Pressure Difference',
                           'from 10 Minutes Ago']

    legend_dps = np.linspace( MAX_PRESSURE_HPA_DIFF, -MAX_PRESSURE_HPA_DIFF, 5 )
    for legend_dp in legend_dps:
        legend_marker = mlines.Line2D( [], [],
                                       color = color_map( legend_dp ),
                                       marker = 'o',
                                       linestyle = 'None',
                                       markersize = 5,
                                       label = f'{legend_dp:+4.2f} hPa' )
        legend_items.append( legend_marker )

    legend_text_lines = ['Visualized by ',
                         'Takashi Nakamoto']
    for legend_text_line in legend_text_lines:
        legend_text = mlines.Line2D( [], [], marker = 'None', linestyle = 'None',
                                     label = legend_text_line)
        legend_items.append( legend_text )

    ax.legend( handles = legend_items,
               title = '\n'.join(legend_title_lines),
               loc = 'center right' )

    date_time_str = date_time.strftime('%Y%m%d_%H%M')
    fig_output_filename = FIG_OUTPUT_DIR / f'{date_time_str}.png'
    fig.savefig( fig_output_filename )
    print( f'Generated {fig_output_filename}' )

    for p in points:
        pp = p.pop(0)
        pp.remove()
        del pp

#    for l in lines.collections:
#        l.remove()
#        del l
