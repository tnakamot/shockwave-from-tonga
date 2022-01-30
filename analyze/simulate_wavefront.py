#!/usr/bin/env python3

#
# Program to draw the wavefront of the shockwave from Hunga Tonga.
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
from multiprocessing import Pool
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from geopy.distance import geodesic
from tqdm import tqdm

from .common import *
from .wavefront import *

def draw_frame(
        minutes_from_eruption,
        width_inch,
        height_inch,
        dpi,
        wavefront_lines,
        projection_method,
        central_longitude,
        central_latitude,
        dump_dir,
        show_message = True,
):
    time_since_eruption = timedelta( minutes = minutes_from_eruption )
    date_time = ERUPTION_TIME + time_since_eruption

    fig = plt.figure( figsize = ( width_inch, height_inch ), dpi = dpi )
    if projection_method == 'PlateCarree':
        projections = [ ccrs.PlateCarree( central_longitude = central_longitude ) ]
    elif projection_method == 'Mollweide':
        projections = [ ccrs.Mollweide( central_longitude = central_longitude ) ]
    elif projection_method == 'Robinson':
        projections = [ ccrs.Robinson( central_longitude = central_longitude ) ]
    elif projection_method == 'Orthographic':
        projections = [ ccrs.Orthographic( central_longitude = central_longitude,
                                           central_latitude  = central_latitude ),
                        ccrs.Orthographic( central_longitude = 180 + central_longitude,
                                           central_latitude  = - central_latitude )]

    for projection_i, projection in enumerate( projections ):
        ax = fig.add_subplot( 1, len( projections ), projection_i + 1, projection = projection )

        ax.set_global()
        ax.stock_img()
        ax.add_feature( cfea.COASTLINE )
        ax.add_feature( cfea.OCEAN )
        ax.add_feature( cfea.LAND )
    
        legend_items = []

        for wavefront_line in wavefront_lines:
            distance_m = wavefront_line.travel_speed_m_s * ( date_time - ERUPTION_TIME ).total_seconds()
            lines = draw_wavefront( ax,
                                    distance       = geodesic( meters = distance_m ),
                                    projection     = projection,
                                    wavefront_line = wavefront_line )

            legend_wavefront_line = mlines.Line2D( [], [] )
            legend_wavefront_line.update_from( lines[0][0] )
            legend_wavefront_line.set_label( f'Speed = {wavefront_line.travel_speed_m_s:d} m/s' )
            legend_items.append( legend_wavefront_line )

        # Draw legend
        if projection_i == 0:
            ax.legend( handles = legend_items,
                       loc = 'upper left',
                       fontsize = 'x-small' )
        
    # Set title
    hours_since_eruption   = int( time_since_eruption.total_seconds() / 3600 )
    minutes_since_eruption = int( ( time_since_eruption.total_seconds() % 3600 ) / 60 )

    title  = 'Estimated wavefront from Hunga Tonga\n'
    title += f'{hours_since_eruption:3d}:{minutes_since_eruption:02d} since the eruption '
    title += '[' + date_time.astimezone( timezone.utc ).strftime('%Y-%m-%d %H:%M:%S (UTC)') + ']'
    fig.suptitle( title )

    img = fig2img( fig, pad_inches = 0.1 )
    plt.close( fig )
    
    if dump_dir is not None:
        img_output_filename = dump_dir / f'simulated_wavefront_{minutes_from_eruption:05d}_minutes.png'
        img.save( img_output_filename )

        # tqdm.write doesn't work as expected if it is called in multiprocessing environment
        if show_message:
            tqdm.write( f'Generated {img_output_filename}' )

    return np.array( img )

def generate_animation(
        output_dir,
        wavefront_lines,
        interval_minutes,
        duration_minutes,
        width_inch,
        height_inch,
        dpi,
        animation_speed_seconds_hours,
        projection_method,
        central_longitude,
        central_latitude,
        multi_process = True,
        dump_dir = None,
):
    output_dir.mkdir( parents = True, exist_ok = True )
    if dump_dir:
        dump_dir.mkdir( parents = True, exist_ok = True )
    
    animation_data = AnimationData()
    simulation_range = range( 0, duration_minutes + interval_minutes, interval_minutes )
    seconds_per_frame = animation_speed_seconds_hours * interval_minutes / 60.0

    if multi_process:
        pool = Pool( processes = len( os.sched_getaffinity(0) ) )

        one_process = partial(
            draw_frame,
            width_inch        = width_inch,
            height_inch       = height_inch,
            dpi               = dpi,
            wavefront_lines   = wavefront_lines,
            projection_method = projection_method,
            central_longitude = central_longitude,
            central_latitude  = central_latitude,
            dump_dir          = dump_dir,
            show_message      = False,
        )

        img_arrays = list( tqdm( pool.imap( one_process, list( simulation_range ) ),
                                 total = len( simulation_range ) ) )

        for img_array in img_arrays:
            animation_data.append(
                Image.fromarray( img_array ),
                duration_ms = int( seconds_per_frame * 1000.0 ),
            )
    else:
        for minutes_from_eruption in tqdm( simulation_range ):
            img_array = draw_frame(
                minutes_from_eruption,
                width_inch,
                height_inch,
                dpi,
                wavefront_lines,
                projection_method,
                central_longitude,
                central_latitude,
                dump_dir,
                show_message = True,
            )
        
            animation_data.append(
                Image.fromarray( img_array ),
                duration_ms = int( seconds_per_frame * 1000.0 )
            )
       
    animation_data.add_cover( text = 'Estimated wavefront of the shockwave\nfrom Hunga Tonga',
                              fontsize = 24,
                              duration_ms = 2000 )
    gif_output_filename = output_dir / f'wavefront_simulation.gif'
    animation_data.save_gif( gif_output_filename )
    print( f'Generated {gif_output_filename}' )

    mp4_output_filename = output_dir / f'wavefront_simulation.mp4'
    animation_data.save_mp4( mp4_output_filename )
    print( f'Generated {mp4_output_filename}' )

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description = 'This program generates animation of estimated wavefront from Hunga Tonga.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        '--outdir',
        dest = 'outdir',
        default = DEFAULT_OUTPUT_DIR,
        help = 'Output directory.',
    )
    parser.add_argument(
        '--projection',
        dest    = 'projection_method',
        choices = ['PlateCarree', 'Mollweide', 'Robinson', 'Orthographic'],
        default = 'PlateCarree',
        help    = 'Map projection method.',
    )
    parser.add_argument(
        '--longitude',
        dest    = 'central_longitude',
        type    = float,
        default = 180.0,
        limit   = ( -180.0, 180.0 ),
        action  = RangeArgument,
        help    = 'Central longitude of the map.',
    )
    parser.add_argument(
        '--latitude',
        dest    = 'central_latitude',
        type    = float,
        default = 0.0,
        limit   = ( -90.0, 90.0 ),
        action  = RangeArgument,
        help    = 'Central latitude of the map. This parameter is effective only when the projection method is Orthographic.',
    )
    parser.add_argument(
        '--interval',
        dest    = 'interval_minutes',
        type    = int,
        default = 10,
        limit   = ( 1, 60 ),
        action  = RangeArgument,
        help    = 'Simulation interval in minutes.',
    )
    parser.add_argument(
        '--duration',
        dest    = 'duration_hours',
        type    = int,
        default = 216,
        limit   = ( 1, 216 ),
        action  = RangeArgument,
        help    = 'Simulation duration in hours.',
    )
    parser.add_argument(
        '--animation-speed',
        dest    = 'animation_speed',
        type    = float,
        default = 1.0,
        limit   = ( 0.1, 10.0 ),
        action  = RangeArgument,
        help    = 'Animation speed in seconds/hours.',
    )
    parser.add_argument(
        '--width',
        dest    = 'width_inch',
        type    = float,
        default = 10.0,
        limit   = ( 1.0, 100.0 ),
        action  = RangeArgument,
        help    = 'Width of the map in inches.',
    )
    parser.add_argument(
        '--height',
        dest    = 'height_inch',
        type    = float,
        default = 5.0,
        limit   = ( 1.0, 100.0 ),
        action  = RangeArgument,
        help    = 'Height of the map in inches.',
    )
    parser.add_argument(
        '--dpi',
        dest    = 'dpi',
        type    = float,
        default = 96.0,
        limit   = ( 1.0, 1000.0 ),
        action  = RangeArgument,
        help    = 'Dot per inch.',
    )
    parser.add_argument(
        '--dump',
        dest = 'dumpdir',
        metavar = 'DUMPDIR',
        help = 'Dump intermediate results to the specified directory.',
    )
    parser.add_argument(
        '--single-process',
        action = 'store_true',
        help    = 'Use a single process instead of multi processes. Useful for debugging.',
    )
    
    return parser

def main():
    '''Main function'''

    parser = create_argument_parser()
    args = parser.parse_args()

    generate_animation(
        output_dir       = args.outdir,
        wavefront_lines  = [
            WavefrontLine( travel_speed_m_s = 320, color = '#000000' ),
            WavefrontLine( travel_speed_m_s = 315, color = '#FF00BF' ),
            WavefrontLine( travel_speed_m_s = 310, color = '#8000FF' ),
            WavefrontLine( travel_speed_m_s = 305, color = '#FFBF00' ),
            WavefrontLine( travel_speed_m_s = 300, color = '#FF0000' ),
        ],
        interval_minutes  = args.interval_minutes,
        duration_minutes  = args.duration_hours * 60,
        width_inch        = args.width_inch,
        height_inch       = args.height_inch,
        dpi               = args.dpi,
        animation_speed_seconds_hours = args.animation_speed,
        projection_method = args.projection_method,
        central_longitude = args.central_longitude,
        central_latitude  = args.central_latitude,
        multi_process     = not args.single_process,
        dump_dir          = Path( args.dumpdir ) if args.dumpdir else None,
    )
   
if __name__ == "__main__":
    main()
