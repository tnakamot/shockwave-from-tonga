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

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from geopy.point import Point
from tqdm import tqdm

from .datasource.asos1min import BarometricPressureDataSourceAsosOneMinute
from .datasource.jma import BarometricPressureDataSourceJMA
from .database import \
    AnalysisDatabase, \
    BarometricPressureMonitoringStation, \
    BarometricPressureMonitoringStationData
from .common import *

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
    'asos1min': ( BarometricPressureDataSourceAsosOneMinute(),
                  MapRange( Point( latitude = 25.0, longitude = - 65.0 ),
                            Point( latitude = 50.0, longitude = -125.0 ) ) ),
    'jma'     : ( BarometricPressureDataSourceJMA(),
                  MapRange( Point( latitude = 24.0, longitude = 120.0 ),
                            Point( latitude = 46.0, longitude = 160.0 ) ) ),
}

def generate_animation(
        database,
        data_source,
        map_range,
        shockwave_i,
        output_dir,
        dump_dir = None
):
    output_dir.mkdir( parents = True, exist_ok = True )
    if dump_dir:
        dump_dir.mkdir( parents = True, exist_ok = True )

    start_time, end_time = map_range.get_time_range( shockwave_i )
    projection = ccrs.PlateCarree()
    animation_data = AnimationData( dump_mode = (dump_dir is not None) )

    # TODO: adjust accordingly
    max_pressure_diff_hPa_minute = 0.05
    pressure_diff_color_map = lambda dp: cm.rainbow( dp / max_pressure_diff_hPa_minute / 2.0 + 0.5 )
    markersize = 3
        
    records = []
    recorded_timestamps = set()
    stations = BarometricPressureMonitoringStation.get_all_stations_in_datasource( database, data_source )
    tqdm.write( f'{len(stations)} stations were found' )
   
    for station in stations:
        data = BarometricPressureMonitoringStationData( database, station, start_time, end_time )
        pressure_diff_hPa_minute, timestamp = data.pressure_diff_hPa_minute_with_timestamp( interpolate = None )
        for i in range( len( timestamp ) ):
            records.append( (station, timestamp[i], pressure_diff_hPa_minute[i] ) )
            recorded_timestamps.add( timestamp[i] )

    tqdm.write( f'{len(records)} records were found' )

    # TODO: parameterize the figure size
    FIG_SIZE = (10, 5)
    fig = plt.figure( figsize = FIG_SIZE )
    recorded_timestamps = sorted( list( recorded_timestamps ) )
    for timestamp in tqdm( recorded_timestamps ):
        ax = fig.add_subplot( 1, 1, 1, projection = projection )
        ax.set_extent( ( map_range.point1.longitude, map_range.point2.longitude,
                         map_range.point1.latitude,  map_range.point2.latitude  ) )
        ax.add_feature( cfea.OCEAN, color = '#00FFFF' )
        ax.add_feature( cfea.LAND,  color = '#32CD32' )
        
        matched_records = [ record for record in records if record[1] == timestamp ]
        for station, _t, pressure_diff_hPa_minute in matched_records:
            ax.plot( station.longitude_deg,
                     station.latitude_deg,
                     'o',
                     transform = projection,
                     markersize = markersize,
                     color = pressure_diff_color_map( pressure_diff_hPa_minute ) )

        # TODO: draw estimated wavefront

        # Generate legend.
        legend_items = []
        utc_timestamp = timestamp.astimezone( timezone.utc )
        legend_title_lines = [ f'{utc_timestamp.strftime("%Y-%m-%d %H:%M")} (UTC)',
                               'Barometric Pressure',
                               'Time Derivative']

        legend_dps = np.linspace( max_pressure_diff_hPa_minute, -max_pressure_diff_hPa_minute, 5 )
        for legend_dp in legend_dps:
            legend_marker = mlines.Line2D( [], [],
                                           color = pressure_diff_color_map( legend_dp ),
                                           marker = 'o',
                                           linestyle = 'None',
                                           markersize = markersize,
                                           label = f'{legend_dp:+4.2f} hPa/minute' )
            legend_items.append( legend_marker )

        legend_text_lines = ['Image by ',
                             'Takashi Nakamoto']

        legend = ax.legend(
            handles = legend_items,
            title = '\n'.join(legend_title_lines),
            loc = 'center right',
            fontsize = 'x-small'
        )
        legend.get_title().set_fontsize('x-small')
        
        # Add this frame to the animation.
        if dump_dir is not None:
            timestamp_str = timestamp.strftime('%Y%m%d_%H%M')
            fig_output_filename = dump_dir / f'{data_source.name}_shockwave_{timestamp_str}.png'

        animation_data.append( fig,
                               duration_ms = 200, #TODO: parameterize duration
                               dump_path = dump_dir )

        if dump_dir is not None:
            tqdm.write( f'Generated {fig_output_filename}' )
            
        fig.clear()

    animation_data.add_cover(
        text = f'{ordinal( shockwave_i + 1 )} shockwave from Hunga Tonga',
        fontsize = 24,
        duration_ms = 2000
    )
    gif_output_filename = output_dir / f'{data_source.name}_shockwave_{shockwave_i}.gif'
    animation_data.save_gif( gif_output_filename )
    tqdm.write( f'Generated {gif_output_filename}' )
    
    mp4_output_filename = output_dir / f'{data_source.name}_shockwave_{shockwave_i}.mp4'
    animation_data.save_mp4( mp4_output_filename )
    tqdm.write( f'Generated {mp4_output_filename}' )

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description = 'This program generates animation of barometric pressure difference on a local map.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    # TODO: this can be moved to common.py.
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
        'source',
        choices = DATA_SOURCES.keys(),
        help = 'Data source name.',
    )
    
    return parser

def main():
    '''Main function'''

    parser = create_argument_parser()
    args = parser.parse_args()

    database = AnalysisDatabase( dbfile_path = args.db )
    data_source, map_range = DATA_SOURCES[ args.source ]

    generate_animation(
        database      = database,
        data_source   = data_source,
        map_range     = map_range,
        shockwave_i   = 0,
        output_dir    = args.outdir, 
        dump_dir      = Path( args.dumpdir ) if args.dumpdir else None
    )

if __name__ == "__main__":
    main()
