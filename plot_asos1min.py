#!/usr/bin/env python3

#
# Program to visualize the shockwave from Tonga using data provided by
# ASOS1MIN network.
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

import os
from datetime import datetime, timedelta, timezone
from functools import partial
from math import isnan, nan
from multiprocessing import Pool
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from geopy.distance import geodesic
from tqdm import tqdm

from common import *

# TODO: move this class out to common.py
class PressureData:
    def __init__(self, latitude_deg, longitude_deg, distance_km, date_time, pressure_hPa, pressure_hPa_diff):
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.distance_km   = distance_km
        self.date_time     = date_time
        self.pressure_hPa  = pressure_hPa
        self.pressure_hPa_diff = pressure_hPa_diff

def read_asos1min_pressure_data(filename):
    f = open( filename, 'r', encoding = 'utf-8' )

    lines = f.readlines()
    latitude_deg  = float( lines[4].split(',')[1] )
    longitude_deg = float( lines[5].split(',')[1] )
    distance_km   = geodesic( HUNGA_TONGA_COORD, ( latitude_deg, longitude_deg ) ).km

    records = []
    pressure_hPa_diffs = []
    
    for line in lines[9:]:
        fields = [ r.strip() for r in line.split(',') ]
        date_time = datetime.fromisoformat( fields[0] )
        pressure_inchHg   = float( fields[1] )
        pressure_hPa      = pressure_inchHg * 33.863886666667

        if records:
            time_diff = date_time - records[-1].date_time
            time_diff_minutes = time_diff.total_seconds() / 60
            pressure_hPa_diff = ( pressure_hPa - records[-1].pressure_hPa ) / time_diff_minutes
        else:
            pressure_hPa_diff = nan

        records.append( PressureData( latitude_deg, longitude_deg, distance_km, date_time, pressure_hPa, pressure_hPa_diff ) )
        pressure_hPa_diffs.append( pressure_hPa_diff )

    f.close()

    return records

def draw_usa_map(fig, projection):
    START_LATITUDE_DEG  =  25.0
    END_LATITUDE_DEG    =  50.0
    START_LONGITUDE_DEG = - 65.0
    END_LONGITUDE_DEG   = -125.0

    ax = fig.add_subplot( 1, 1, 1, projection = projection )
    ax.set_extent( ( START_LONGITUDE_DEG, END_LONGITUDE_DEG,
                     START_LATITUDE_DEG, END_LATITUDE_DEG ),
                   projection )
    ax.add_feature( cfea.OCEAN, color = OCEAN_COLOR )
    ax.add_feature( cfea.LAND,  color = LAND_COLOR )
    ax.add_feature( cfea.STATES )
    return ax

def generate_animation(fig,
                       gif_output_filepath,
                       start_time,
                       end_time,
                       pressure_hPa_diff_max,
                       pressure_hPa_diff_min,
                       title,
                       records):
    DUMP_OUTPUT_DIR = Path( 'figure_asos1min/dump' )
    DUMP_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
    
    max_pressure_hPa_diff = max( abs( pressure_hPa_diff_max ), abs( pressure_hPa_diff_min ) )
    pressure_diff_color_map = lambda dp: cm.rainbow( dp / max_pressure_hPa_diff / 2.0 + 0.5 )

    dump_mode = True
    animation_data = AnimationData( dump_mode = dump_mode )

    recorded_date_times = sorted( list( set( [ record.date_time.astimezone( timezone.utc )
                                               for record in records ] ) ) )
    date_times = [ date_time for date_time in recorded_date_times if start_time <= date_time <= end_time ]
    projection = ccrs.PlateCarree()
    
    for date_time in tqdm( date_times[::10] ): # TODO generate every one minute
        ax = draw_usa_map( fig, projection )

        matched_records = [ record for record in records
                            if record.date_time == date_time and not isnan( record.pressure_hPa_diff ) ]
        
        # Plot the pressure difference data in the map.
        for matched_record in matched_records:
            ax.plot( matched_record.longitude_deg,
                     matched_record.latitude_deg,
                     'o',
                     transform = projection,
                     markersize = 5,
                     color = pressure_diff_color_map( matched_record.pressure_hPa_diff ) )

        # Generate legend.
        legend_items = []
        legend_title_lines = [ f'{date_time.strftime("%Y-%m-%d %H:%M")} (UTC)',
                               'Pressure Difference',
                               'from 1 minute ago' ]

        legend_dps = np.linspace( max_pressure_hPa_diff, -max_pressure_hPa_diff, 5 )
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
                   loc = 'lower right' )
            
        date_time_str = date_time.strftime('%Y%m%d_%H%M')
        fig_output_filename = DUMP_OUTPUT_DIR / f'{date_time_str}.png'

        animation_data.append( fig,
                               duration_ms = 5,
                               dump_path = fig_output_filename )

        if dump_mode:
            tqdm.write( f'Generated {fig_output_filename}' )

        fig.clear()

    animation_data.add_cover( text = title,
                              fontsize = 24,
                              duration_ms = 2000 )
    animation_data.save_gif( gif_output_filepath )
    tqdm.write( f'Generated {gif_output_filepath}' )

def generate_time_distance_scatter_plot(output_filename,
                                        start_time,
                                        end_time,
                                        pressure_hPa_diff_max,
                                        pressure_hPa_diff_min,
                                        title,
                                        records):
    distances_km = []
    hours_since_eruption = []
    pressure_hPa_diffs = []

    for record in tqdm( records ):
        if isnan( record.pressure_hPa_diff ):
            continue
        if record.date_time < start_time or record.date_time > end_time:
            continue
        if record.distance_km < 8000:
            continue
        
        distances_km.append( record.distance_km )
        hours_since_eruption.append( ( record.date_time - ERUPTION_TIME ).total_seconds() / 3600 )
        pressure_hPa_diffs.append( record.pressure_hPa_diff )

    min_hours = np.min( hours_since_eruption )
    max_hours = np.max( hours_since_eruption )
    min_distance_km = np.min( distances_km )
    max_distance_km = np.max( distances_km )
    
    fig = plt.figure()
    ax = fig.add_subplot( 1, 1, 1 )
    im = ax.scatter( hours_since_eruption,
                     distances_km,
                     c = pressure_hPa_diffs,
                     cmap = cm.rainbow,
                     linewidth = 0,
                     vmin = pressure_hPa_diff_min,
                     vmax = pressure_hPa_diff_max,
                     s = 1 )
    ax.set_xlabel( 'Time since eruption [hours]' )
    ax.set_xlim( min_hours, max_hours )
    ax.set_ylabel( 'Distance from Hunga Tonga [km]' )
    ax.set_ylim( min_distance_km, max_distance_km )
    yticks = ax.get_yticks()
    yticks = [min_distance_km] + list( yticks ) + [max_distance_km]
    ax.set_yticks( yticks )
    ax.set_ylim( min_distance_km, max_distance_km )
    ax.grid( True, which = 'major',
             linewidth = 1, linestyle = '--', color = '#808080' )
    fig.colorbar( im, ax = ax,
                  label = 'Pressure difference from 1 minute ago [hPa]' )
    ax.set_title( title )

    fig.savefig( output_filename )

def estimate_la_arrival_times():
    LA_COORD = Point( latitude  = 34.05,
                      longitude = -118.25 )
    HUNGA_TONGA_TO_LA = geodesic( HUNGA_TONGA_COORD, LA_COORD )

    # Travel distance of the 1st, 2nd, 3rd, 4th and 5th shockwaves.
    travel_distances = [
        HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE - HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE + HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 2 - HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 2 + HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 3 - HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 3 + HUNGA_TONGA_TO_LA,
    ]

    # Estimated travel speed of shockwave eastwards.
    TRAVEL_SPEED_M_S_TO_EAST = 320 # [m/s]

    # Estaimted travel speed of shockwave westwards.
    TRAVEL_SPEED_M_S_TO_WEST = 307 # [m/s]

    # Estimated arrival time of the shockwaves at Log Angeles.
    estimated_la_arrival_times = []
    for shockwave_i, travel_distance in enumerate( travel_distances ):
        travel_speed_m_s = TRAVEL_SPEED_M_S_TO_EAST if shockwave_i % 2 == 0 else TRAVEL_SPEED_M_S_TO_WEST
        estimated_la_arrival_times.append(
            ERUPTION_TIME + timedelta( seconds = travel_distance.m / travel_speed_m_s )
        )
    return estimated_la_arrival_times
    
def main():
    # Read pressure data from CSV files.
    print('Reading CSV files...')
    SKIP_ANIMATION_GENERATION = False
    
    DATA_INPUT_DIR = Path( 'data_asos1min' )
    csv_files = list( DATA_INPUT_DIR.glob( '*.csv' ) )
    records = [ record for csv_file in tqdm( csv_files ) for record in read_asos1min_pressure_data( csv_file ) ]

    FIG_OUTPUT_DIR = Path( 'figure_asos1min' )
    FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

    # Generate GIF animation.
    FIG_US_MAP_SIZE = (14.4, 6)
    fig_us_map = plt.figure( figsize = FIG_US_MAP_SIZE )

    start_time = datetime( 2022, 1, 15, 11, 30, tzinfo = timezone.utc )
    end_time   = datetime( 2022, 1, 15, 18, 00, tzinfo = timezone.utc )

    for shockwave_i, estimated_la_arrival_time in enumerate( tqdm( estimate_la_arrival_times() ) ):
        if shockwave_i % 2 == 0:
            start_time = estimated_la_arrival_time - timedelta( hours = 0.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 7.5 )
        else:
            start_time = estimated_la_arrival_time - timedelta( hours = 5.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 2.5 )

        if shockwave_i < 2:
            pressure_hPa_diff_max =  0.25
            pressure_hPa_diff_min = -0.25
        elif shockwave_i < 4:
            pressure_hPa_diff_max =  0.15
            pressure_hPa_diff_min = -0.15
        else:
            pressure_hPa_diff_max =  0.05
            pressure_hPa_diff_min = -0.05

        time_distance_filepath = FIG_OUTPUT_DIR / f'time_distance_shockwave_{shockwave_i}.png'
        title  = f'{ordinal(shockwave_i+1)} shockwave from Hunga Tonga\n'
        title +=  'US ASOS one minute interval pressure data'
        generate_time_distance_scatter_plot( output_filename = time_distance_filepath,
                                             start_time      = start_time,
                                             end_time        = end_time,
                                             pressure_hPa_diff_max = pressure_hPa_diff_max,
                                             pressure_hPa_diff_min = pressure_hPa_diff_min,
                                             title           = title,
                                             records         = records )
        tqdm.write( f'Generated {time_distance_filepath}' )

        if SKIP_ANIMATION_GENERATION:
            continue

        # Generating animation on the map takes very long time.
        gif_animation_filepath = FIG_OUTPUT_DIR / f'us_map_with_shockwave_{shockwave_i}.gif'
        title = f'{ordinal(shockwave_i+1)} shockwave from Hunga Tonga'
        generate_animation( fig                   = fig_us_map,
                            gif_output_filepath   = gif_animation_filepath,
                            start_time            = start_time,
                            end_time              = end_time,
                            pressure_hPa_diff_max = pressure_hPa_diff_max,
                            pressure_hPa_diff_min = pressure_hPa_diff_min,
                            title                 = title,
                            records               = records )

if __name__ == "__main__":
    main()
