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

import os
import sys
sys.path.append( os.path.join( os.path.dirname( __file__ ) ) )

import io
from datetime import datetime, timedelta, timezone
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from geopy.distance import geodesic
from tqdm import tqdm

from common import *

def draw_wavefront(fig, minutes_from_eruption, travel_speed_m_s):
    time_since_eruption = timedelta( minutes = minutes_from_eruption )
    date_time = ERUPTION_TIME + time_since_eruption
    bearings = np.linspace(-180, 180, 360)
    distance = geodesic( meters = travel_speed_m_s * ( date_time - ERUPTION_TIME ).total_seconds() )
    wavefront_points = [ distance.destination( point = HUNGA_TONGA_COORD, bearing = b ) for b in bearings ]

    wavefront_latitude_deg  = [ p[0] for p in wavefront_points ]
    wavefront_longitude_deg = [ p[1] for p in wavefront_points ]

    ax = add_world_map( fig, ccrs.PlateCarree() )

    div_js = np.where( np.abs( np.diff( wavefront_longitude_deg ) ) > 180 )[0] + 1
    div_js = np.append( div_js, len( wavefront_longitude_deg ) )
    div_j_start = 0
    for div_j_end in div_js:
        ax.plot( wavefront_longitude_deg[div_j_start:div_j_end],
                 wavefront_latitude_deg[div_j_start:div_j_end],
                 transform = ccrs.PlateCarree(),
                 color = 'black' )
        div_j_start = div_j_end

    hours_since_eruption   = int( time_since_eruption.total_seconds() / 3600 )
    minutes_since_eruption = int( ( time_since_eruption.total_seconds() % 3600 ) / 60 )
        
    title  = 'Estimated wavefront from Hunga Tonga\n'
    title += f'{hours_since_eruption:3d}:{minutes_since_eruption:02d} since the eruption'
    title += '[' + date_time.astimezone( timezone.utc ).strftime('%Y-%m-%d %H:%M:%S (UTC)') + ']\n'
    title += f'(Assumed travel speed = {travel_speed_m_s:.1f} m/s)'
    ax.set_title( title )


def main():
    # Estimated travel speed of shock wave to draw.
    TRAVEL_SPEED_M_S = 310 # [m/s]

    OUTPUT_DIR = Path( 'figure_wavefront_simulation' )
    OUTPUT_DIR.mkdir( parents = True, exist_ok = True )
    
    SIMULATION_INTERVAL_MINUTES = 10
    SIMULATION_DURATION_MINUTES = 72 * 60

    fig = plt.figure( figsize = (10, 5) )
    dump_mode = True
    animation_data = AnimationData( dump_mode = dump_mode )
    simulation_range = range( 0,
                              SIMULATION_DURATION_MINUTES + SIMULATION_INTERVAL_MINUTES,
                              SIMULATION_INTERVAL_MINUTES )

    for minutes_from_eruption in tqdm( simulation_range ):
        draw_wavefront( fig, minutes_from_eruption, TRAVEL_SPEED_M_S )
        
        dump_path = OUTPUT_DIR / f'wavefrom_simulation_{minutes_from_eruption:05d}.png'
        
        animation_data.append( fig,
                               duration_ms = 100,
                               dump_path = dump_path )
        if dump_mode:
            tqdm.write( f'Generated {dump_path}' )

        fig.clear()
        
        
    animation_data.add_cover( text = 'Estimated wavefront of the shockwave\nfrom Hunga Tonga',
                              fontsize = 24,
                              duration_ms = 2000 )
    gif_output_filename = OUTPUT_DIR / f'wavefront_simulation.gif'
    animation_data.save_gif( gif_output_filename )
    print( f'Generated {gif_output_filename}' )

if __name__ == "__main__":
    main()
