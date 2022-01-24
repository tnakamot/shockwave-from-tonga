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
                       pressure_diff_max_hPa_minute,
                       pressure_diff_min_hPa_minute,
                       title,
                       records):
    DUMP_OUTPUT_DIR = Path( 'figure_asos1min/dump' )
    DUMP_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
    
    max_pressure_hPa_diff = max( abs( pressure_diff_max_hPa_minute ), abs( pressure_diff_min_hPa_minute ) )
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


class AsosStations:
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

    def ids_in_range( self,
                      min_distance_km = 0,
                      max_distance_km = 40000 ):
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


def generate_time_distance_scatter_plot(sqlite3_connection,
                                        shockwave_i,
                                        start_time,
                                        end_time,
                                        pressure_diff_max_hPa_minute,
                                        pressure_diff_min_hPa_minute):
    MIN_DISTANCE_KM = 8000

    stations = AsosStations( sqlite3_connection )
    
    hours_since_eruption = []
    distance_km = []
    pressure_diff_hPa_minute = []

    filtered_hours_since_eruption = []
    filtered_distance_km = []
    filtered_pressure_diff_hPa_minute = []

    distance_km_for_line = []
    peak_hour_for_line = []

    cursor = sqlite3_connection.cursor()

    for station_id in stations.ids_in_range( min_distance_km = MIN_DISTANCE_KM ):
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
    ''', (station_id,
          start_time.astimezone( timezone.utc ),
          end_time.astimezone( timezone.utc ) ) )

        minutes = []
        pressure_hPas = []

        previous_timestamp = None
        previous_pressure_hPa = None
        for row in cursor.fetchall():
            current_timestamp    = datetime.fromisoformat( row[0] )
            current_pressure_hPa = row[1]

            seconds_since_eruption = ( current_timestamp - ERUPTION_TIME ).total_seconds()
            minutes.append( seconds_since_eruption / 60 )
            pressure_hPas.append( current_pressure_hPa )

            if previous_timestamp and \
               ( ( current_timestamp - previous_timestamp ) == timedelta( minutes = 1 ) ) :
                hours_since_eruption.append( seconds_since_eruption / 3600 )
                pressure_diff_hPa_minute.append( current_pressure_hPa - previous_pressure_hPa )
                distance_km.append( stations.distance_km( station_id ) )
                
            previous_timestamp    = current_timestamp
            previous_pressure_hPa = current_pressure_hPa


        if not minutes:
            continue
        
        start_minutes  = ( start_time - ERUPTION_TIME ).total_seconds() / 60
        end_minutes    = ( end_time - ERUPTION_TIME ).total_seconds() / 60
        interp_minutes = np.arange( start_minutes, end_minutes, 1 )
        interp_pressure_hPas = np.interp( interp_minutes, minutes, pressure_hPas )

        # High-pass filter
        #
        # filter_weight = firwin( numtaps   = 151,
        #                         cutoff    = [0.03, 0.2],
        #                         window    = 'hamming',
        #                         pass_zero = 'bandpass',
        #                         fs        = 1.0 )
        # filtered_pressure_hPas = np.convolve( interp_pressure_hPas, filter_weight, 'same' )
        # filtered_pressure_diff_hPa_minute.extend( np.diff( filtered_pressure_hPas ) )
        # filtered_hours_since_eruption.extend( interp_minutes[1:] / 60 )
        # filtered_distance_km.extend( [ stations.distance_km( station_id ) ] * ( len(interp_minutes) - 1 ) )

        # Hilbert transform
        # interp_pressure_diff_hPa_minute = np.diff( interp_pressure_hPas )
        # filtered_pressure_diff_hPa_minute.extend( np.abs( hilbert( interp_pressure_diff_hPa_minute ) ) )
        # filtered_hours_since_eruption.extend( interp_minutes[1:] / 60 )
        # filtered_distance_km.extend( [ stations.distance_km( station_id ) ] * ( len(interp_minutes) - 1 ) )

        # Hilbert transform & peak detection
        interp_pressure_diff_hPa_minute = np.diff( interp_pressure_hPas )
        envelope = np.abs( hilbert( interp_pressure_diff_hPa_minute ) )
        peak = np.zeros( len( envelope ) )
        peak[ np.argmax( envelope ) ] = 1.0
        filtered_pressure_diff_hPa_minute.extend( peak )
        filtered_hours_since_eruption.extend( interp_minutes[1:] / 60 )
        filtered_distance_km.extend( [ stations.distance_km( station_id ) ] * ( len(interp_minutes) - 1 ) )

        distance_km_for_line.append( stations.distance_km( station_id ) )
        peak_hour_for_line.append( interp_minutes[ np.argmax( envelope ) ] / 60 )

    p_fitted = np.polyfit( distance_km_for_line, peak_hour_for_line, 1 )
    peak_hour_for_fitted_line = np.poly1d( p_fitted )
    estimated_travel_speed_m_s = 1.0 / p_fitted[0] * 1000.0 / 3600.0

    min_hours = np.min( hours_since_eruption )
    max_hours = np.max( hours_since_eruption )
    min_distance_km = np.min( distance_km )
    max_distance_km = np.max( distance_km )

    images = {}
    fig = plt.figure()
    for filtered in [False, True]:
        ax = fig.add_subplot( 1, 1, 1 )
        x = filtered_hours_since_eruption if filtered else hours_since_eruption
        y = filtered_distance_km if filtered else distance_km
        z = filtered_pressure_diff_hPa_minute if filtered else pressure_diff_hPa_minute
        vmin = pressure_diff_min_hPa_minute if filtered else pressure_diff_min_hPa_minute
        vmax = pressure_diff_max_hPa_minute if filtered else pressure_diff_max_hPa_minute
        im = ax.scatter( x, y, c = z,
                         cmap = cm.rainbow,
                         linewidth = 0,
                         vmin = vmin,
                         vmax = vmax,
                         s = 1 )
        if filtered or not filtered:
            ax.plot( peak_hour_for_fitted_line( distance_km_for_line ),
                     distance_km_for_line,
                     color = 'black', linestyle = '-', linewidth = 1 )
            ax.text( 0.5, 0.5,
                     f'Estimated travel speed\n{estimated_travel_speed_m_s:.1f} m/s',
                     horizontalalignment = 'left',
                     verticalalignment = 'center',
                     transform = ax.transAxes )

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

        title  = f'{ordinal(shockwave_i+1)} shockwave from Hunga Tonga\n'
        title +=  'US ASOS one minute interval pressure data'
        if filtered:
            title += ' (max of hilbert abs)'
        ax.set_title( title )

        FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
        if filtered:
            filepath = FIG_OUTPUT_DIR / f'filtered_time_distance_shockwave_{shockwave_i}.png'
        else:
            filepath = FIG_OUTPUT_DIR / f'time_distance_shockwave_{shockwave_i}.png'
            
        fig.savefig( filepath, bbox_inches = 'tight', pad_inches = 0.05 )
        tqdm.write( f'Generated {filepath}.' )

        images[ 'filtered' if filtered else 'raw' ] = fig2img( fig )
        fig.clear()
        
    return images

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
    def __init__(self,
                 shockwave_i,
                 start_time,
                 end_time,
                 pressure_diff_max_hPa_minute,
                 pressure_diff_min_hPa_minute):
        self.shockwave_i = shockwave_i
        self.start_time  = start_time
        self.end_time    = end_time
        self.pressure_diff_max_hPa_minute = pressure_diff_max_hPa_minute
        self.pressure_diff_min_hPa_minute = pressure_diff_min_hPa_minute

def generate_shockwave_parameters(shockwave_num):
    estimated_la_arrival_times = estimate_la_arrival_times()
    params = []

    for shockwave_i in range(shockwave_num):
        if shockwave_i < 2:
            pressure_diff_max_hPa_minute =  0.25
            pressure_diff_min_hPa_minute = -0.25
        elif shockwave_i < 4:
            pressure_diff_max_hPa_minute =  0.15
            pressure_diff_min_hPa_minute = -0.15
        else:
            pressure_diff_max_hPa_minute =  0.05
            pressure_diff_min_hPa_minute = -0.05

        if shockwave_i % 2 == 0:
            pressure_diff_max_hPa_minute /= 2
            pressure_diff_min_hPa_minute /= 2

        estimated_la_arrival_time = estimated_la_arrival_times[shockwave_i]
        if shockwave_i % 2 == 0:
            start_time = estimated_la_arrival_time - timedelta( hours = 0.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 7.5 )
        else:
            start_time = estimated_la_arrival_time - timedelta( hours = 4.5 )
            end_time   = estimated_la_arrival_time + timedelta( hours = 3.5 )

        param = ShockwaveParameter( shockwave_i = shockwave_i,
                                    start_time  = start_time,
                                    end_time    = end_time,
                                    pressure_diff_max_hPa_minute = pressure_diff_max_hPa_minute,
                                    pressure_diff_min_hPa_minute = pressure_diff_min_hPa_minute )
        params.append( param )
    return params

def process_one_shockwave( shockwave_param ):

    sqlite3_connection = sqlite3.connect( ASOS1MIN_SQLITE3_DATABASE )

    images = generate_time_distance_scatter_plot( sqlite3_connection = sqlite3_connection,
                                                  shockwave_i        = shockwave_param.shockwave_i,
                                                  start_time         = shockwave_param.start_time,
                                                  end_time           = shockwave_param.end_time,
                                                  pressure_diff_max_hPa_minute = shockwave_param.pressure_diff_max_hPa_minute,
                                                  pressure_diff_min_hPa_minute = shockwave_param.pressure_diff_min_hPa_minute )

    sqlite3_connection.close()
    
    return images

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

    pool = Pool( processes = len( os.sched_getaffinity(0) ) )
    time_distance_images = list( tqdm( pool.imap( process_one_shockwave, shockwave_params ),
                                       total = len( shockwave_params ) ) )

    # Generate combined time vs distance images.
    for image_type, filename in [ ('raw',      'time_distance_shockwaves.png'),
                                  ('filtered', 'filtered_time_distance_shockwaves.png') ]:
        images = [ img_d[image_type] for img_d in time_distance_images ]
        combined_time_distance_image = combine_images( images, padding = 10, portrait = False )
        
        FIG_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
        combined_filepath = FIG_OUTPUT_DIR / filename
        combined_time_distance_image.save( combined_filepath )
        print( f'Generated {combined_filepath}' )

    # Move the code below
        
        # if SKIP_ANIMATION_GENERATION:
        #     continue

        # # Generating animation on the map takes very long time.
        # gif_animation_filepath = FIG_OUTPUT_DIR / f'us_map_with_shockwave_{shockwave_i}.gif'
        # title = f'{ordinal(shockwave_i+1)} shockwave from Hunga Tonga'
        # generate_animation( fig                   = fig_us_map,
        #                     gif_output_filepath   = gif_animation_filepath,
        #                     start_time            = start_time,
        #                     end_time              = end_time,
        #                     pressure_diff_max_hPa_minute = pressure_diff_max_hPa_minute,
        #                     pressure_diff_min_hPa_minute = pressure_diff_min_hPa_minute,
        #                     title                 = title,
        #                     records               = records )

if __name__ == "__main__":
    main()
