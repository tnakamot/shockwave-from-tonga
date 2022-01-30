#!/usr/bin/env python3

#
# Program to visualize the shockwave from Tonga using data.
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

import argparse
from datetime import datetime, timedelta, timezone
from functools import partial

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from geopy.point import Point
from tqdm import tqdm

from .datasource.asos1min import BarometricPressureDataSourceAsosOneMinute
from .datasource.jma import BarometricPressureDataSourceJMA
from .database import \
    AnalysisDatabase, \
    BarometricPressureMonitoringStation, \
    BarometricPressureMonitoringStationData
from .common import *
from .wavefront import *

class MapRange:
    # Roughly estimated propagation speed of shockwave.
    ESTIMATED_PROPAGATION_SPEED_M_S = 315 # [m/s]
    
    def __init__( self, point1, point2 ):
        self.point1 = point1
        self.point2 = point2

    def get_time_range( self, shockwave_i ):
        point3 = Point( self.point1.latitude, self.point2.longitude )
        point4 = Point( self.point2.latitude, self.point1.longitude )
        points = [ self.point1, self.point2, point3, point4 ]

        time_since_eruption = []
        for point in points:
            if shockwave_i % 2 == 0:
                origin_coord = HUNGA_TONGA_COORD
            else:
                origin_coord = ANTIPODE_COORD

            distance = EARTH_CIRCUMFERENCE * shockwave_i / 2.0 + geodesic( origin_coord, point )
            time_since_eruption.append( timedelta( seconds = distance.m / self.ESTIMATED_PROPAGATION_SPEED_M_S ) )

        START_TIME_MARGIN = timedelta( hours = 1 )
        END_TIME_MARGIN   = timedelta( hours = 3 )

        return ( ERUPTION_TIME + min( time_since_eruption ) - START_TIME_MARGIN,
                 ERUPTION_TIME + max( time_since_eruption ) + END_TIME_MARGIN )

DATA_SOURCES = {
    'asos1min': {
        'data_source': BarometricPressureDataSourceAsosOneMinute(),
        'map_range'  : MapRange( Point( latitude = 23.0, longitude = - 65.0 ),
                                 Point( latitude = 48.0, longitude = -125.0 ) ),
        'acknowledgement': \
'''Data from US ASOS sites
through Iowa Environmental Mesonet
of Iowa State University''',
        'default_width_inch' : 14.4,
        'default_height_inch':  6.0,
        'default_dpi'        : 96.0,
        'legend_loc'         : 'lower right',
    },
    'jma'     : {
        'data_source': BarometricPressureDataSourceJMA(),
        'map_range'  : MapRange( Point( latitude = 23.0, longitude = 120.0 ),
                                 Point( latitude = 45.0, longitude = 160.0 ) ),
        'acknowledgement': 'Data from Japan Meteorological Agency',
        'default_width_inch' : 10.0,
        'default_height_inch':  5.0,
        'default_dpi'        : 96.0,
        'legend_loc'         : 'center right',
    },
}

def generate_snapshot(
        records,
        data_source,
        map_range,
        width_inch,
        height_inch,
        dpi,
        legend_loc,
        dump_dir,
        show_wavefront,
        acknowledgement,
        max_pressure_diff_hPa_minute,
        show_message = True,
):
    timestamp = records[0][1]
    
    projection = ccrs.PlateCarree()
    markersize = 4
    pressure_diff_color_map = lambda dp: cm.seismic( dp / max_pressure_diff_hPa_minute / 2.0 + 0.5 )
    
    fig = plt.figure( figsize = ( width_inch, height_inch ), dpi = dpi )
    ax = fig.add_subplot( 1, 1, 1, projection = projection )
    ax.set_extent( ( map_range.point1.longitude, map_range.point2.longitude,
                     map_range.point1.latitude,  map_range.point2.latitude  ) )
    ax.add_feature( cfea.OCEAN, color = '#FFFFFF' )
    ax.add_feature( cfea.LAND,  color = '#C0C0C0' )
    ax.add_feature( cfea.BORDERS )
    ax.add_feature( cfea.STATES )

    for station, _t, pressure_diff_hPa_minute in records:
        ax.plot( station.longitude_deg,
                 station.latitude_deg,
                 'o',
                 transform = projection,
                 markersize = markersize,
                 markeredgewidth = 0.5,
                 markeredgecolor = 'black',
                 color = pressure_diff_color_map( pressure_diff_hPa_minute ) )

    # Draw estimated wavefront.
    wavefront_lines = [
        WavefrontLine( travel_speed_m_s = 320, color = '#000000' ),
        WavefrontLine( travel_speed_m_s = 315, color = '#FF00BF' ),
        WavefrontLine( travel_speed_m_s = 310, color = '#8000FF' ),
        WavefrontLine( travel_speed_m_s = 305, color = '#FFBF00' ),
        WavefrontLine( travel_speed_m_s = 300, color = '#FF0000' ),
    ]

    legend_wavefront_lines = []
    for wavefront_line in wavefront_lines:
        if not show_wavefront:
            break
                
        distance_m = wavefront_line.travel_speed_m_s * ( timestamp - ERUPTION_TIME ).total_seconds()
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
    utc_timestamp = timestamp.astimezone( timezone.utc )
    legend_title_lines = [ f'{utc_timestamp.strftime("%Y-%m-%d %H:%M")} (UTC)',
                            'Barometric Pressure Time Derivative']
       
    legend_dps = np.linspace( max_pressure_diff_hPa_minute, -max_pressure_diff_hPa_minute, 5 )
    for legend_dp in legend_dps:
        legend_marker = mlines.Line2D( [], [],
                                       color = pressure_diff_color_map( legend_dp ),
                                       marker = 'o',
                                       linestyle = 'None',
                                       markersize = markersize,
                                       markeredgewidth = 0.5,
                                       markeredgecolor = 'black',
                                       label = f'{legend_dp:+5.3f} hPa/minute' )
        legend_items.append( legend_marker )

    for legend_wavefront_line in legend_wavefront_lines:
        legend_items.append( legend_wavefront_line )
            
    legend_text_lines  = ['Image by Takashi Nakamoto']
    legend_text_lines += acknowledgement.split('\n')
    for legend_text_line in legend_text_lines:
        legend_text = mlines.Line2D( [], [], marker = 'None', linestyle = 'None',
                                         label = legend_text_line)
        legend_items.append( legend_text )

    legend = ax.legend(
        handles = legend_items,
        title = '\n'.join(legend_title_lines),
        loc = legend_loc,
        fontsize = 'x-small'
    )
    legend.get_title().set_fontsize('x-small')

    img = fig2img( fig )
    plt.close( fig )

    if dump_dir is not None:
        timestamp_str = timestamp.strftime('%Y%m%d_%H%M')
        img_output_filename = dump_dir / f'{data_source.name}_shockwave_{timestamp_str}.png'
        img.save( img_output_filename )

        # tqdm.write doesn't work as expected if it is called in multiprocessing environment
        if show_message:
            tqdm.write( f'Generated {img_output_filename}' )
    
    return np.array( img )

def generate_animation(
        database,
        data_source,
        map_range,
        shockwave_i,
        output_dir,
        width_inch,
        height_inch,
        dpi,
        max_pressure_diff_hPa_minute,
        legend_loc = 'center right',
        dump_dir = None,
        show_wavefront = False,
        acknowledgement = '',
        multi_process = True,
):
    # Output directory parameters
    output_dir.mkdir( parents = True, exist_ok = True )
    if dump_dir:
        dump_dir.mkdir( parents = True, exist_ok = True )

    # Animation parameters
    start_time, end_time = map_range.get_time_range( shockwave_i )
    animation_data = AnimationData()

    # Obtain stations in the range
    stations = BarometricPressureMonitoringStation.get_all_stations_in_datasource( database, data_source )
    tqdm.write( f'{len(stations)} stations were found' )

    # Obtain pressure time derivative data for each station
    records = []
    recorded_timestamps = set()
    for station in stations:
        data = BarometricPressureMonitoringStationData( database, station, start_time, end_time )
        pressure_diff_hPa_minute, timestamp = data.pressure_diff_hPa_minute_with_timestamp( interpolate = None )
        for i in range( len( timestamp ) ):
            records.append( (station, timestamp[i], pressure_diff_hPa_minute[i] ) )
            recorded_timestamps.add( timestamp[i] )

    tqdm.write( f'{len(records)} records were found' )

    recorded_timestamps = sorted( list( recorded_timestamps ) )
    
    dt_minutes = (recorded_timestamps[1] - recorded_timestamps[0]).total_seconds() / 60
    duration_ms = 40 * dt_minutes

    if multi_process:
        pool = Pool( processes = len( os.sched_getaffinity(0) ) )

        one_process = partial(
            generate_snapshot,
            data_source = data_source,
            map_range   = map_range,
            width_inch  = width_inch,
            height_inch = height_inch,
            dpi         = dpi,
            legend_loc  = legend_loc,
            dump_dir    = dump_dir,
            show_wavefront  = show_wavefront,
            acknowledgement = acknowledgement,
            max_pressure_diff_hPa_minute = max_pressure_diff_hPa_minute,
            show_message    = False,
        )

        split_records = []
        for timestamp in recorded_timestamps:
            split_records.append( [ record for record in records if record[1] == timestamp ] )
            
        img_arrays = list( tqdm( pool.imap( one_process, split_records ),
                                 total = len( split_records ) ) )
        for img_array in img_arrays:
            animation_data.append(
                Image.fromarray( img_array ),
                duration_ms = duration_ms,
            )
        
    else:
        for timestamp in tqdm( recorded_timestamps ):
            img_array = generate_snapshot(
                [ record for record in records if record[1] == timestamp ],
                data_source,
                map_range,
                width_inch,
                height_inch,
                dpi,
                legend_loc,
                dump_dir,
                show_wavefront,
                acknowledgement,
                max_pressure_diff_hPa_minute = max_pressure_diff_hPa_minute,
            )
        
            animation_data.append(
                Image.fromarray( img_array ),
                duration_ms = duration_ms,
            )

    animation_data.add_cover(
        text = f'{ordinal( shockwave_i + 1 )} shockwave from Hunga Tonga',
        fontsize = 24,
        duration_ms = 2000
    )
    wavefront_suffix = '_with_wavefront' if show_wavefront else ''
    gif_output_filename = output_dir / f'{data_source.name}_shockwave_{shockwave_i + 1}{wavefront_suffix}.gif'
    animation_data.save_gif( gif_output_filename )
    tqdm.write( f'Generated {gif_output_filename}' )
    
    mp4_output_filename = output_dir / f'{data_source.name}_shockwave_{shockwave_i + 1}{wavefront_suffix}.mp4'
    animation_data.save_mp4( mp4_output_filename )
    tqdm.write( f'Generated {mp4_output_filename}' )

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description = 'This program generates animation of barometric pressure difference on a local map.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    # TODO: "--db" and "source" this can be moved to common.py.
    parser.add_argument(
        '--db',
        default = AnalysisDatabase.DEFAULT_DBFILE_PATH,
        help = 'Path to the SQLite3 database file.',
    )
    parser.add_argument(
        '--outdir',
        dest = 'outdir',
        default = DEFAULT_OUTPUT_DIR,
        help = 'Output directory.',
    )
    parser.add_argument(
        '--dump',
        dest = 'dumpdir',
        metavar = 'DUMPDIR',
        help = 'Dump intermediate results to the specified directory.',
    )
    parser.add_argument(
        '--wavefront',
        dest   = 'show_wavefront',
        action = 'store_true',
        help   = 'Switch to draw estimated wavefront.',
    )
    parser.add_argument(
        '--width',
        dest   = 'width_inch',
        type   = float,
        help   = 'Width of the map in inches.',
    )
    parser.add_argument(
        '--height',
        dest   = 'height_inch',
        type   = float,
        help   = 'Height of the map in inches.',
    )
    parser.add_argument(
        '--dpi',
        dest   = 'dpi',
        type   = float,
        help   = 'Dot per inch.',
    )
    parser.add_argument(
        '--scale',
        type    = float,
        default = 0.1,
        help    = 'Maximum scale of pressure time derivative in hPa/minute.',
    )
    parser.add_argument(
        '--single_process',
        action = 'store_true',
        help    = 'Use a single process instead of multi processes. Useful for debugging.',
    )
    parser.add_argument(
        'source',
        choices = DATA_SOURCES.keys(),
        help = 'Data source name.',
    )
    parser.add_argument(
        'shockwave_i',
        choices = range(1,9),
        type    = int,
        help    = 'Generate animation of i-th shockwave.',
    )
    
    return parser

def main():
    '''Main function'''

    parser = create_argument_parser()
    args = parser.parse_args()

    database = AnalysisDatabase( dbfile_path = args.db )
    width_inch = args.width_inch if args.width_inch else DATA_SOURCES[ args.source ][ 'default_width_inch' ]
    height_inch = args.height_inch if args.height_inch else DATA_SOURCES[ args.source ][ 'default_height_inch' ]
    dpi = args.dpi if args.dpi else DATA_SOURCES[ args.source ][ 'default_dpi' ]

    generate_animation(
        database        = database,
        data_source     = DATA_SOURCES[ args.source ][ 'data_source' ],
        map_range       = DATA_SOURCES[ args.source ][ 'map_range' ],
        shockwave_i     = args.shockwave_i - 1,
        output_dir      = args.outdir,
        width_inch      = width_inch,
        height_inch     = height_inch,
        dpi             = dpi,
        max_pressure_diff_hPa_minute = args.scale,
        multi_process   = not args.single_process,
        legend_loc      = DATA_SOURCES[ args.source ][ 'legend_loc' ],
        dump_dir        = Path( args.dumpdir ) if args.dumpdir else None,
        show_wavefront  = args.show_wavefront,
        acknowledgement = DATA_SOURCES[ args.source ][ 'acknowledgement' ],
    )

if __name__ == "__main__":
    main()
