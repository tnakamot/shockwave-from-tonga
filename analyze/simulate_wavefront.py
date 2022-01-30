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
        width_inch,
        height_inch,
        dpi,
        minutes_from_eruption,
        wavefront_lines,
        dump_dir,
        show_message = True,
):
    fig = plt.figure( figsize = ( width_inch, height_inch ), dpi = dpi )
    projection = ccrs.PlateCarree()
    
    ax = fig.add_subplot( 1, 1, 1, projection = projection )
    ax.set_extent( ( -180, 180, -90, 90 ), projection )

    # TODO: use better colors
    ax.add_feature( cfea.OCEAN, color = '#8080FF' )
    ax.add_feature( cfea.LAND,  color = '#80FF80' )
    
    time_since_eruption = timedelta( minutes = minutes_from_eruption )
    date_time = ERUPTION_TIME + time_since_eruption
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
    ax.legend( handles = legend_items,
               loc = 'upper left',
               fontsize = 'x-small' )
        
    # Set title
    hours_since_eruption   = int( time_since_eruption.total_seconds() / 3600 )
    minutes_since_eruption = int( ( time_since_eruption.total_seconds() % 3600 ) / 60 )

    title  = 'Estimated wavefront from Hunga Tonga\n'
    title += f'{hours_since_eruption:3d}:{minutes_since_eruption:02d} since the eruption '
    title += '[' + date_time.astimezone( timezone.utc ).strftime('%Y-%m-%d %H:%M:%S (UTC)') + ']'
    ax.set_title( title )

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
        multi_process = True,
        dump_dir = None,
):
    output_dir.mkdir( parents = True, exist_ok = True )
    if dump_dir:
        dump_dir.mkdir( parents = True, exist_ok = True )
    
    animation_data = AnimationData()
    simulation_range = range( 0, duration_minutes + interval_minutes, interval_minutes )
    seconds_per_frame = animation_speed_seconds_hours * interval_minutes / 60.0

    for minutes_from_eruption in tqdm( simulation_range ):
        img_array = draw_frame(
            width_inch,
            height_inch,
            dpi,
            minutes_from_eruption,
            wavefront_lines,
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
        '--dump',
        dest = 'dumpdir',
        metavar = 'DUMPDIR',
        help = 'Dump intermediate results to the specified directory.',
    )
    parser.add_argument(
        '--width',
        dest    = 'width_inch',
        type    = float,
        default = 10.0,
        help    = 'Width of the map in inches.',
    )
    parser.add_argument(
        '--height',
        dest    = 'height_inch',
        type    = float,
        default = 5.0,
        help    = 'Height of the map in inches.',
    )
    parser.add_argument(
        '--dpi',
        dest    = 'dpi',
        type    = float,
        default = 96.0,
        help    = 'Dot per inch.',
    )
    parser.add_argument(
        '--interval',
        dest    = 'interval_minutes',
        type    = int,
        default = 10,
        help    = 'Simulation interval in minutes.',
    )
    parser.add_argument(
        '--duration',
        dest    = 'duration_hours',
        type    = int,
        default = 216,
        help    = 'Simulation duration in hours.',
    )
    parser.add_argument(
        '--animation-speed',
        dest    = 'animation_speed',
        type    = float,
        default = 1.0,
        help    = 'Animation speed in seconds/hours.',
    )
    parser.add_argument(
        '--single_process',
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
        interval_minutes = args.interval_minutes,
        duration_minutes = args.duration_hours * 60,
        width_inch       = args.width_inch,
        height_inch      = args.height_inch,
        dpi              = args.dpi,
        animation_speed_seconds_hours = args.animation_speed,
        multi_process    = not args.single_process,
        dump_dir         = Path( args.dumpdir ) if args.dumpdir else None,
    )
   
if __name__ == "__main__":
    main()
