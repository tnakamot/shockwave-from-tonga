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

import itertools
import os
import sqlite3
from datetime import datetime, timedelta, timezone
from functools import partial
from math import ceil
from multiprocessing import Pool
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from geopy.distance import geodesic
from scipy.signal import firwin, hilbert
from scipy.interpolate import interp1d
from tqdm import tqdm

from common import *

# Change this parameter to False to generate shockwave animation over the map os USA,
# but it takes a very long time.
SKIP_ANIMATION_GENERATION = True 

# The size of USA map.
FIG_US_MAP_SIZE = (14.4, 6)

# Output directory.
FIG_OUTPUT_DIR = Path( 'figure_asos1min' )

# Table name in the SQLite3 database.
STATION_TABLE_NAME = 'station'
PRESSURE_DATA_TABLE_NAME = 'pressure_data'

# Minimum distance in km from Hunga Tonga for time distance scatter plot
MIN_DISTANCE_KM = 8000

class AsosStations:
    '''This class handles database access to the ASOS station information.'''
    
    def __init__(self, sqlite3_connection):
        self._sqlite3_connection = sqlite3_connection
        self._load()

    def _load(self):
        cursor = self._sqlite3_connection.cursor()
        cursor.execute( f'''
SELECT
    id,
    distance_km
FROM
    {STATION_TABLE_NAME}
ORDER BY
    distance_km
''')
        rows = cursor.fetchall()
        self._station_ids = [ station_id for station_id, distance_km in rows ]
        self._distance_km = { station_id: distance_km for station_id, distance_km in rows }
        
    def all_ids(self):
        '''Get all station IDs of the ASOS1MIN network in ascending order of distance from Hunga Tonga.

        Returns:
        --------
        list
            List of station IDs (e.g., AGC). Typically, it is three or four letter string.
            See the following web page to see what are the station IDs:
            https://mesonet.agron.iastate.edu/request/asos/1min.phtml
        '''
        return self._station_ids

    def ids_in_range( self, min_distance_km = 0, max_distance_km = 40000 ):
        '''Get station IDs of the ASOS1MIN network stations in the specified range of distance from Hunga Tonga.

        Parameters:
        -----------
        min_distance_km: float
            Minimum distance from Hunga Tonga in km.

        max_distance_km: float
            Maximum distance from Hunga Tonga in km.

        Returns:
        --------
        list
            List of station IDs (e.g., AGC). They are in ascending order of distance from Hunga Tonga.
        '''
        return [ station_id for station_id in self._station_ids
                 if min_distance_km <= self._distance_km[station_id] <= max_distance_km ]

    def distance_km( self, station_id ):
        '''Return the distance of the station from Hunga Tonga.

        Parameters:
        -----------
        station_id: str
            Station ID.

        Returns:
        --------
        float
            Distance from Hunga Tonga in km.
        '''
        return self._distance_km[ station_id ]

class StationPressureData:
    '''This class handles database access to the pressure data of the specified station.'''

    def __init__( self, sqlite3_connection, station_id, start_time, end_time ):
        self._sqlite3_connection = sqlite3_connection
        self._station_id = station_id
        self._start_time = start_time.astimezone( timezone.utc )
        self._end_time   =   end_time.astimezone( timezone.utc )
        self._load()

    def _load( self ):
        cursor = self._sqlite3_connection.cursor()
        cursor.execute( f'''
SELECT
    timestamp,
    pressure_hPa
FROM
    {PRESSURE_DATA_TABLE_NAME}
WHERE
    station_id = ? AND
    timestamp >= ? AND
    timestamp <  ?
ORDER BY
    timestamp ASC
        ''', ( self._station_id,
               self._start_time,
               self._end_time ) )
        rows = cursor.fetchall()
        self._timestamps   = [ datetime.fromisoformat( row[0] ) for row in rows ]
        self._pressure_hPa = [ row[1] for row in rows ]

    def is_empty( self ):
        '''Returns True if there is no pressure data or only one record in the specified date and time range.'''
        return len( self._timestamps ) <= 1

    def _get_pressure_hPa( self, interpolate = None ):
        minutes_since_eruption = np.array( [ ( t - ERUPTION_TIME ).total_seconds() / 60
                                             for t in self._timestamps ] )
        if not interpolate:
            return self._pressure_hPa, minutes_since_eruption

        start_minutes = int( ( self._start_time - ERUPTION_TIME ).total_seconds() / 60 ) + 1
        end_minutes   = int( ( self._end_time   - ERUPTION_TIME ).total_seconds() / 60 )
        interpolated_minutes_since_eruption = np.arange( start_minutes, end_minutes, 1 )
        f_interpolate_pressure_hPa = interp1d( minutes_since_eruption,
                                               self._pressure_hPa,
                                               kind = interpolate,
                                               copy = False,
                                               assume_sorted = True )
        extrapolate_start_n = int(minutes_since_eruption[0] - start_minutes)
        extrapolate_end_n   = int(end_minutes - minutes_since_eruption[-1]) 
        interp_start_i = extrapolate_start_n
        interp_end_i   = len( interpolated_minutes_since_eruption ) - extrapolate_end_n
        
        interpolated_pressure_hPa = np.concatenate( (
            np.full( extrapolate_start_n, self._pressure_hPa[0] ),
            f_interpolate_pressure_hPa( interpolated_minutes_since_eruption[ interp_start_i : interp_end_i ] ),
            np.full( extrapolate_end_n, self._pressure_hPa[-1] )
        ) )
        
        return interpolated_pressure_hPa, interpolated_minutes_since_eruption

    def _pressure_diff_hPa_minute_interp( self, interpolate ):
        pressure_hPa, minutes_since_eruption = self._get_pressure_hPa( interpolate )
        _pressure_diff_hPa_minute = np.diff( pressure_hPa )
        return _pressure_diff_hPa_minute, minutes_since_eruption[1:]
    
    def pressure_diff_hPa_minute( self, interpolate = None ):
        if interpolate:
            return self._pressure_diff_hPa_minute_interp( interpolate )

        previous_timestamp = None
        previous_pressure_hPa = None
        minutes_since_eruption = []
        _pressure_diff_hPa_minute = []
        for i, current_timestamp in enumerate( self._timestamps ):
            current_pressure_hPa = self._pressure_hPa[ i ]

            if previous_timestamp and \
               ( ( current_timestamp - previous_timestamp ) == timedelta( minutes = 1 ) ):
                seconds_since_eruption = ( current_timestamp - ERUPTION_TIME ).total_seconds()
                minutes_since_eruption.append( seconds_since_eruption / 60 )
                _pressure_diff_hPa_minute.append( current_pressure_hPa - previous_pressure_hPa )

            previous_timestamp    = current_timestamp
            previous_pressure_hPa = current_pressure_hPa

        return _pressure_diff_hPa_minute, minutes_since_eruption

def configure_time_distance_scatter_ax(
        ax,
        hours_since_eruption,
        distance_km,
        additional_title,
        shockwave_i,
):

    min_distance_km = min( distance_km )
    max_distance_km = max( distance_km )
    
    ax.set_xlabel( 'Time since eruption [hours]' )
    ax.set_xlim( min( hours_since_eruption), max( hours_since_eruption ) )
    ax.set_ylabel( 'Distance from Hunga Tonga [km]' )
    ax.set_ylim( min_distance_km, max_distance_km )
    
    yticks = ax.get_yticks()
    yticks = [ min_distance_km ] + list( yticks ) + [ max_distance_km ]
    ax.set_yticks( yticks )
    ax.set_ylim( min_distance_km, max_distance_km )
    
    ax.grid( True, which = 'major',
             linewidth = 1, linestyle = '--', color = '#808080' )

    title  = f'{ordinal(shockwave_i+1)} shockwave from Hunga Tonga\n'
    title +=  'Based on US ASOS one minute interval pressure data\n'
    title +=  additional_title
    ax.set_title( title )
    
def generate_time_distance_scatter_plot(
        stations,
        data,
        shockwave_i,
        scale_mode,
        interpolate = None
):
    hours_since_eruption = []
    pressure_diff_hPa_minute = []
    distance_km = []

    for station_id in stations.ids_in_range( min_distance_km = MIN_DISTANCE_KM ):
        if data[station_id].is_empty():
            continue

        ps, ms = data[station_id].pressure_diff_hPa_minute( interpolate )
        hours_since_eruption.extend( [ m / 60.0 for m in ms ] )
        pressure_diff_hPa_minute.extend( ps )
        distance_km.extend( [ stations.distance_km( station_id ) ] * len( ps ) )

    max_pressure_diff_hPa_minute = 0.05
    if scale_mode == 'compare':
        if shockwave_i == 0 or shockwave_i == 1:
            max_pressure_diff_hPa_minute = 0.20
        elif shockwave_i == 2 or shockwave_i == 3:
            max_pressure_diff_hPa_minute = 0.10
    elif scale_mode == 'best' and 0 <= shockwave_i <= 3:
        max_pressure_diff_hPa_minute = (0.2, 0.25, 0.075, 0.15)[shockwave_i]

    fig = plt.figure()
    ax = fig.add_subplot( 1, 1, 1 )
    im = ax.scatter(
        hours_since_eruption,
        distance_km,
        c = pressure_diff_hPa_minute,
        cmap = cm.rainbow,
        linewidth = 0,
        vmin = - max_pressure_diff_hPa_minute,
        vmax = + max_pressure_diff_hPa_minute,
        s = 1
    )

    configure_time_distance_scatter_ax( ax,
                                        hours_since_eruption,
                                        distance_km,
                                        'Pressure difference from 1 minute ago',
                                        shockwave_i )
    fig.colorbar( im, ax = ax,
                  label = 'Pressure difference from 1 minute ago [hPa]' )

    img = fig2img( fig )
    fig.clear()
        
    return img

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
        EARTH_CIRCUMFERENCE * 4 - HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 4 + HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 5 - HUNGA_TONGA_TO_LA,
        EARTH_CIRCUMFERENCE * 5 + HUNGA_TONGA_TO_LA,
    ]

    # Estimated travel speed of shockwave eastwards.
    TRAVEL_SPEED_M_S_TO_EAST = 320 # [m/s]

    # Estaimted travel speed of shockwave westwards.
    TRAVEL_SPEED_M_S_TO_WEST = 310 # [m/s]

    # Estimated arrival time of the shockwaves at Log Angeles.
    estimated_la_arrival_times = []
    for shockwave_i, travel_distance in enumerate( travel_distances ):
        travel_speed_m_s = TRAVEL_SPEED_M_S_TO_EAST if shockwave_i % 2 == 0 else TRAVEL_SPEED_M_S_TO_WEST
        estimated_la_arrival_time = ERUPTION_TIME + timedelta( seconds = travel_distance.m / travel_speed_m_s )
        estimated_la_arrival_times.append( estimated_la_arrival_time )
    return estimated_la_arrival_times

class ShockwaveParameter:
    def __init__( self,
                  shockwave_i,
                  start_time,
                  end_time
    ):
        self.shockwave_i = shockwave_i
        self.start_time  = start_time
        self.end_time    = end_time

def generate_shockwave_parameters(shockwave_num):
    estimated_la_arrival_times = estimate_la_arrival_times()
    params = []

    for shockwave_i in range(shockwave_num):
        estimated_la_arrival_time = estimated_la_arrival_times[shockwave_i]
        if shockwave_i % 2 == 0:
            start_time = estimated_la_arrival_time - timedelta( hours = 0.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 7.5 )
        else:
            start_time = estimated_la_arrival_time - timedelta( hours = 4.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 3.5 )

        param = ShockwaveParameter(
            shockwave_i = shockwave_i,
            start_time  = start_time,
            end_time    = end_time,
        )
        
        params.append( param )
    return params

class VisualizedShockwave:
    def __init__( self, shockwave_i, chart_type, image ):
        self.shockwave_i = shockwave_i
        self.chart_type  = chart_type
        self.image       = image

def visualize_one_shockwave( params ):
    shockwave_param, chart_type = params

    sqlite3_connection = sqlite3.connect( ASOS1MIN_SQLITE3_DATABASE )
    stations = AsosStations( sqlite3_connection )
    data = { station_id: \
             StationPressureData( sqlite3_connection,
                                  station_id,
                                  shockwave_param.start_time,
                                  shockwave_param.end_time )
             for station_id in stations.all_ids() }

    chart_generators = {
        'raw_same_scale': lambda: generate_time_distance_scatter_plot(
            stations    = stations,
            data        = data,
            shockwave_i = shockwave_param.shockwave_i,
            scale_mode  = 'same',
            interpolate = None,
        ),
        'raw_best_scale': lambda: generate_time_distance_scatter_plot(
            stations    = stations,
            data        = data,
            shockwave_i = shockwave_param.shockwave_i,
            scale_mode  = 'best',
            interpolate = None,
        ),
        'raw_comparison': lambda: generate_time_distance_scatter_plot(
            stations    = stations,
            data        = data,
            shockwave_i = shockwave_param.shockwave_i,
            scale_mode  = 'compare',
            interpolate = None,
        ),
        'interpolated': lambda: generate_time_distance_scatter_plot(
            stations    = stations,
            data        = data,
            shockwave_i = shockwave_param.shockwave_i,
            scale_mode  = 'compare',
            interpolate = 'linear',
        ),
    }

    image = chart_generators[chart_type]()

    sqlite3_connection.close()
    
    return VisualizedShockwave(
        shockwave_i = shockwave_param.shockwave_i,
        chart_type  = chart_type,
        image       = image,
    )

def combine_images(images, padding = 10, portrait = True):
    one_width  = max( image.width  for image in images )
    one_height = min( image.height for image in images )

    if portrait:
        combined_width  = one_width * 2 + padding * 3
        combined_height = ( one_height + padding ) * ceil( len(images) / 2) + padding
    else:
        combined_width  = ( one_width + padding ) * ceil( len(images) / 2) + padding
        combined_height = one_height * 2 + padding * 3

    combined_image  = Image.new( 'RGB', ( combined_width, combined_height ), 'white' )
    for image_i in range( len(images) ):
        if portrait:
            x = ( one_width + padding )  * ( image_i % 2 ) + padding
            y = ( one_height + padding ) * int( image_i / 2 ) + padding
        else:
            x = ( one_width + padding )  * int( image_i / 2 ) + padding
            y = ( one_height + padding ) * ( image_i % 2 ) + padding
        
        combined_image.paste( images[ image_i ], (x, y) )

    return combined_image

def main():
    shockwave_nums   = 8
    shockwave_params = generate_shockwave_parameters( shockwave_nums )
    chart_type_and_filenames = [
        ('raw_same_scale', 'time_distance_chart_raw_same_scale.png'),
        ('raw_best_scale', 'time_distance_chart_raw_best_scale.png'),
        ('raw_comparison', 'time_distance_chart_raw_scale_for_comparison.png'),
        ('interpolated'  , 'time_distance_chart_interpolated.png'),
    ]

    chart_types_to_generate = [ chart_type for chart_type, filename in chart_type_and_filenames ]
    chart_types_to_generate = ['interpolated']

    params = list( itertools.product( shockwave_params, chart_types_to_generate ) )
    multi_process = True # Turn this switch to False for debugging in a single process mode.
    if multi_process:
        pool = Pool( processes = len( os.sched_getaffinity(0) ) )
        visualized_shockwaves = list( tqdm( pool.imap( visualize_one_shockwave, params ),
                                            total = len( params ) ) )
    else:
        visualized_shockwaves = [ visualize_one_shockwave(param) for param in tqdm( params ) ]

    # Generate combined time vs distance images.
    for chart_type, filename in chart_type_and_filenames:
        if chart_type not in chart_types_to_generate:
            continue
            
        images = [ vs.image for vs in visualized_shockwaves if vs.chart_type == chart_type ]
        combined_image = combine_images( images, padding = 10, portrait = False )
        
        FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
        combined_filepath = FIG_OUTPUT_DIR / filename
        combined_image.save( combined_filepath )
        print( f'Generated {combined_filepath}' )

if __name__ == "__main__":
    main()
