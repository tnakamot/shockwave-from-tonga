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

import io
from datetime import datetime, timedelta, timezone
from math import isnan, nan
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from geopy.distance import geodesic
from geopy.point import Point
from tqdm import tqdm

from common import *

class PressureData:
    def __init__(self, latitude_deg, longitude_deg, date_time, pressure_hPa, pressure_hPa_diff):
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.date_time     = date_time
        self.pressure_hPa  = pressure_hPa
        self.pressure_hPa_diff = pressure_hPa_diff

def read_jma_pressure_data(filename):
    f = open( filename, 'r', encoding = 'utf-8' )
    
    lines = f.readlines()
    latitude_deg  = float( lines[4].split(',')[1] )
    longitude_deg = float( lines[5].split(',')[1] )

    data = []
    pressure_hPa_prev = nan

    for line in lines[9:]:
        fields = [r.strip() for r in line.split(',')]
        date_time    = datetime.fromisoformat( fields[0] )
        pressure_hPa = float( fields[1] )
        pressure_hPa_diff = pressure_hPa - pressure_hPa_prev
        pressure_hPa_prev = pressure_hPa

        data.append( PressureData( latitude_deg, longitude_deg, date_time, pressure_hPa, pressure_hPa_diff ) )
    
    f.close()

    return data

def draw_japan_map(fig, projection):
    START_LATITUDE_DEG  =  24.0
    END_LATITUDE_DEG    =  46.0
    START_LONGITUDE_DEG = 120.0
    END_LONGITUDE_DEG   = 160.0

    ax = fig.add_subplot( 1, 1, 1, projection = projection )
    ax.set_extent( ( START_LONGITUDE_DEG, END_LONGITUDE_DEG,
                     START_LATITUDE_DEG, END_LATITUDE_DEG ),
                   projection )
    ax.add_feature( cfea.OCEAN, color = OCEAN_COLOR )
    ax.add_feature( cfea.LAND,  color = LAND_COLOR )
    return ax


def generate_animation(fig, shockwave_i, start_time, end_time, records):
    FIG_OUTPUT_DIR = Path( 'figure_jma' )
    FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

    wavefront_lines = [
        WavefrontLine( travel_speed_m_s = 320, color = '#000000' ),
        WavefrontLine( travel_speed_m_s = 315, color = '#FF00BF' ),
        WavefrontLine( travel_speed_m_s = 310, color = '#8000FF' ),
        WavefrontLine( travel_speed_m_s = 305, color = '#FFBF00' ),
        WavefrontLine( travel_speed_m_s = 300, color = '#FF0000' ),
    ]

    
    MAX_PRESSURE_HPA_DIFF = 0.5
    pressure_diff_color_map = lambda dp: cm.rainbow( dp / MAX_PRESSURE_HPA_DIFF / 2.0 + 0.5 )

    dump_mode = True
    animation_data = AnimationData( dump_mode = dump_mode )
    
    recorded_date_times = sorted( list( set( [ record.date_time for record in records ] ) ) )
    date_times = [ date_time for date_time in recorded_date_times if start_time <= date_time <= end_time ]
    projection = ccrs.PlateCarree()

    for date_time in tqdm( date_times ):
        ax = draw_japan_map( fig, projection )
        
        matched_records = [ record for record in records if record.date_time == date_time ]

        # Plot the pressure difference data in the map.
        for matched_record in matched_records:
        
            if isnan( matched_record.pressure_hPa_diff ):
                continue

            ax.plot( matched_record.longitude_deg,
                     matched_record.latitude_deg,
                     'o',
                     transform = projection,
                     markersize = 5,
                     color = pressure_diff_color_map( matched_record.pressure_hPa_diff ) )

        # Draw estimated wavefront.
        legend_wavefront_lines = []
        for wavefront_line in wavefront_lines:
            distance_m = wavefront_line.travel_speed_m_s * ( date_time - ERUPTION_TIME ).total_seconds()
            lines = draw_wavefront( ax,
                                    distance       = geodesic( meters = distance_m ),
                                    projection     = projection,
                                    wavefront_line = wavefront_line )
            
            legend_wavefront_line = mlines.Line2D( [], [] )
            legend_wavefront_line.update_from( lines[0][0] )
            legend_wavefront_line.set_label( f'Estimated wavefront ({wavefront_line.travel_speed_m_s:d} m/s)' )
            legend_wavefront_lines.append( legend_wavefront_line )

        # Generate legend.
        legend_items = []
        legend_title_lines = [f'{date_time.strftime("%Y-%m-%d %H:%M")} (JST)',
                              'Barometric Pressure Difference',
                              'from 10 Minutes Ago']

        for legend_wavefront_line in legend_wavefront_lines:
            legend_items.append( legend_wavefront_line )

        legend_dps = np.linspace( MAX_PRESSURE_HPA_DIFF, -MAX_PRESSURE_HPA_DIFF, 5 )
        for legend_dp in legend_dps:
            legend_marker = mlines.Line2D( [], [],
                                           color = pressure_diff_color_map( legend_dp ),
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
                   loc = 'center right',
                   fontsize = 'x-small' )

        date_time_str = date_time.strftime('%Y%m%d_%H%M')
        fig_output_filename = FIG_OUTPUT_DIR / f'{date_time_str}.png'
        
        animation_data.append( fig,
                               duration_ms = 200,
                               dump_path = fig_output_filename )

        if dump_mode:
            tqdm.write( f'Generated {fig_output_filename}')

        fig.clear()

    animation_data.add_cover( text = f'{ordinal( shockwave_i + 1 )} shockwave from Hunga Tonga',
                              fontsize = 24,
                              duration_ms = 2000 )
    gif_output_filename = FIG_OUTPUT_DIR / f'shockwave_{shockwave_i}.gif'
    animation_data.save_gif( gif_output_filename )
    tqdm.write( f'Generated {gif_output_filename}' )

def estimate_kyoto_arrival_times():
    KYOTO_COORD = Point( latitude  = 35.011611,
                         longitude = 135.768111 )
    HUNGA_TONGA_TO_KYOTO = geodesic( HUNGA_TONGA_COORD, KYOTO_COORD )

    # Roughly estimated travel speed of shockwave.
    TRAVEL_SPEED_M_S = 310 # [m/s]

    # Travel distance of the 1st, 2nd, 3rd, 4th and 5th shockwaves.
    travel_distances = [
        HUNGA_TONGA_TO_KYOTO,
        EARTH_CIRCUMFERENCE - HUNGA_TONGA_TO_KYOTO,
        EARTH_CIRCUMFERENCE + HUNGA_TONGA_TO_KYOTO,
        EARTH_CIRCUMFERENCE * 2 - HUNGA_TONGA_TO_KYOTO,
        EARTH_CIRCUMFERENCE * 2 + HUNGA_TONGA_TO_KYOTO,
    ]

    # Estimated arrival time of the shockwaves at Kyoto.
    estimated_kyoto_arrival_times \
        = [ ERUPTION_TIME + timedelta( seconds = d.m / TRAVEL_SPEED_M_S ) for d in travel_distances ]
    return estimated_kyoto_arrival_times
    
def main():
    # Read pressure data from CSV files.
    DATA_INPUT_DIR = Path( 'data_jma' )
    csv_files = DATA_INPUT_DIR.glob( '*.csv' )
    records = [ record for csv_file in csv_files for record in read_jma_pressure_data( csv_file ) ]

    # Generate GIF animation.
    FIG_SIZE = (10, 5)
    fig = plt.figure( figsize = FIG_SIZE )

    for shockwave_i, estimated_kyoto_arrival_time in enumerate( tqdm( estimate_kyoto_arrival_times() ) ):
        start_time = estimated_kyoto_arrival_time - timedelta( hours = 3 )
        end_time   = estimated_kyoto_arrival_time + timedelta( hours = 5 )

        generate_animation( fig, shockwave_i, start_time, end_time, records )

if __name__ == "__main__":
    main()
