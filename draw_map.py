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
from datetime import datetime, timezone, timedelta
from math import nan, isnan
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib.font_manager import findfont, FontProperties
import cartopy.crs as ccrs
import cartopy.feature as cfea
from pathlib import Path
from geopy import distance
from PIL import Image, ImageDraw, ImageFont
import io

TZ_JST = timezone( timedelta( hours = +9 ), name = 'JST' ) # Japan Standard Time
GENERATE_ANIMATION_GIF = True

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
        date_time    = datetime.strptime( fields[0][:-3] + fields[0][-2:], '%Y-%m-%dT%H:%M:%S%z')
        pressure_hPa = float( fields[1] )
        pressure_hPa_diff = pressure_hPa - pressure_hPa_prev
        pressure_hPa_prev = pressure_hPa

        data.append( PressureData( latitude_deg, longitude_deg, date_time, pressure_hPa, pressure_hPa_diff ) )
    
    f.close()

    return data

def fig2img(fig):
    buf = io.BytesIO()
    fig.savefig( buf )
    buf.seek( 0 )
    img = Image.open( buf )
    return img

def ordinal(i):
    if 11 <= (i % 100) <= 13:
        return str(i) + 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd'] + ['th'] * 6
        return str(i) + suffix[i % 10]

# ==============================================================================
# Estimate the arrival time of the shockwave.
# ==============================================================================
hunga_tonga_coord = ( -20.536   , -175.382   )
antipode_coord    = ( - hunga_tonga_coord[0], 180 + hunga_tonga_coord[1] )
kyoto_coord       = (  35.011611,  135.768111)

earth_circumference  = distance.distance( hunga_tonga_coord, antipode_coord ) * 2
hunga_tonga_to_kyoto = distance.distance( hunga_tonga_coord, kyoto_coord )

# Travel distance of the 1st, 2nd, 3rd, 4th and 5th shockwaves.
travel_distances = [
    hunga_tonga_to_kyoto,
    earth_circumference - hunga_tonga_to_kyoto,
    earth_circumference + hunga_tonga_to_kyoto,
    earth_circumference * 2 - hunga_tonga_to_kyoto,
    earth_circumference * 2 + hunga_tonga_to_kyoto,
]

# Hunga Tonga eruption time estimated from the satellite image of Himawari 8
#  https://himawari.asia/himawari8-image.htm?sI=D531106&sClC=ffff00&sTA=true&sTAT=TY&sS=6&sNx=3&sNy=2&sL=-169.171875&sT=-426.8125&wW=1920&wH=969&sD=1642219800000
eruption_time    = datetime( 2022, 1, 15, 13, 10, tzinfo = TZ_JST )

# Estimated travel speed of shockwave.
#  TODO: Justify why it is slower than the speed of sound at 1 atm, 15 deg C (340.5 m/s)
travel_speed_m_s = 310 # [m/s]

# Estimated arrival time of the shockwaves at Kyoto.
estimated_kyoto_arrival_times \
    = [ eruption_time + timedelta( seconds = d.m / travel_speed_m_s ) for d in travel_distances ]

# ==============================================================================
# Read presure data
# ==============================================================================
DATA_INPUT_DIR = Path( 'data_jma' )
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
FIG_OUTPUT_DIR = Path( 'figure_jma' )
FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

MAX_PRESSURE_HPA_DIFF = 0.5
color_map = lambda dp: cm.rainbow( dp / MAX_PRESSURE_HPA_DIFF / 2.0 + 0.5 )

recorded_date_times = sorted( list( set( [ record.date_time for record in records ] ) ) )

for shockwave_i, estimated_kyoto_arrival_time in enumerate( estimated_kyoto_arrival_times ):
    start_time = estimated_kyoto_arrival_time - timedelta( hours = 3 )
    end_time   = estimated_kyoto_arrival_time + timedelta( hours = 5 )
    animation_images = []

    for date_time in recorded_date_times:
        if date_time < start_time or date_time > end_time:
            continue
        
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

        legend_text_lines = ['Image by ',
                             'Takashi Nakamoto']
        for legend_text_line in legend_text_lines:
            legend_text = mlines.Line2D( [], [], marker = 'None', linestyle = 'None',
                                         label = legend_text_line)
            legend_items.append( legend_text )

        ax.legend( handles = legend_items,
                   title = '\n'.join(legend_title_lines),
                   loc = 'center right' )

        if GENERATE_ANIMATION_GIF:
            print(f'Internally generated the map at {date_time.strftime("%Y-%m-%d %H:%M (JST)")}.')
            animation_images.append( fig2img( fig ) )
        else:
            date_time_str = date_time.strftime('%Y%m%d_%H%M')
            fig_output_filename = FIG_OUTPUT_DIR / f'{date_time_str}.png'
            fig.savefig( fig_output_filename, bbox_inches='tight', pad_inches = 0 )
            print( f'Generated {fig_output_filename}' )

        for p in points:
            pp = p.pop(0)
            pp.remove()
            del pp

    if GENERATE_ANIMATION_GIF:
        # Generate the cover page.
        image_size = animation_images[0].size
        cover_image = Image.new( 'RGB', image_size, 'black' )
        cover_draw  = ImageDraw.Draw( cover_image )
        cover_font  = ImageFont.truetype( findfont( FontProperties( family = ['monospace'] ) ),
                                          size = 36 )
        cover_text  = f'{ordinal( shockwave_i + 1 )} shockwave from Hunga Tonga'
        cover_draw.text( ( image_size[0] / 2, image_size[1] / 2 ),
                         cover_text,
                         fill = 'white',
                         font = cover_font,
                         anchor = 'mm',
                         align = 'center' )
        
        # Generate animation GIF.
        duration = [2000] + [200] * len( animation_images )
        gif_output_filename = FIG_OUTPUT_DIR / f'shockwave_{shockwave_i}.gif'
        cover_image.save( gif_output_filename,
                          append_images = animation_images,
                          save_all = True,
                          duration = duration,
                          loop = 0 )
        print( f'Generated {gif_output_filename}' )
