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

import io
from datetime import datetime, timedelta, timezone
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfea
import matplotlib.pyplot as plt
import numpy as np
from geopy.distance import geodesic
from matplotlib.font_manager import FontProperties, findfont
from PIL import Image, ImageDraw, ImageFont

OUTPUT_DIR = Path( 'figure_wavefront_simulation' )
OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

TZ_JST = timezone( timedelta( hours = +9 ), name = 'JST' ) # Japan Standard Time
GENERATE_ANIMATION_GIF = True

OCEAN_COLOR = '#00FFFF'
LAND_COLOR  = '#32CD32'

hunga_tonga_coord   = ( -20.536   , -175.382   )
antipode_coord      = ( - hunga_tonga_coord[0], 180 + hunga_tonga_coord[1] )
earth_circumference = geodesic( hunga_tonga_coord, antipode_coord ) * 2
eruption_time       = datetime( 2022, 1, 15, 13, 10, tzinfo = TZ_JST )
travel_speed_m_s    = 310 # [m/s]

fig = plt.figure( figsize = (10, 5) )

def fig2img(fig):
    buf = io.BytesIO()
    fig.savefig( buf, bbox_inches='tight', pad_inches = 0 )
    buf.seek( 0 )
    img = Image.open( buf )
    return img


animation_images = []

for minutes_from_eruption in range(0, 72 * 60, 10):
    time_since_eruption = timedelta( minutes = minutes_from_eruption )
    date_time = eruption_time + time_since_eruption
    bearings = np.linspace(-180, 180, 360)
    distance = geodesic( meters = travel_speed_m_s * ( date_time - eruption_time ).total_seconds() )
    wavefront_points = [ distance.destination( point = hunga_tonga_coord, bearing = b ) for b in bearings ]

    wavefront_latitude_deg  = [ p[0] for p in wavefront_points ]
    wavefront_longitude_deg = [ p[1] for p in wavefront_points ]

    ax = fig.add_subplot( 1, 1, 1,
                          projection = ccrs.PlateCarree() )
    ax.set_extent( ( -180, 180, -90, 90 ), ccrs.PlateCarree() )
    ax.add_feature( cfea.OCEAN, color = OCEAN_COLOR )
    ax.add_feature( cfea.LAND,  color = LAND_COLOR )

    div_js = np.where( np.abs( np.diff( wavefront_longitude_deg ) ) > 180 )[0] + 1
    div_js = np.append( div_js, len( wavefront_longitude_deg ) )
    div_j_start = 0
    for div_j_end in div_js:
        ax.plot( wavefront_longitude_deg[div_j_start:div_j_end],
                 wavefront_latitude_deg[div_j_start:div_j_end],
                 transform = ccrs.PlateCarree(),
                 color = 'black' )
        div_j_start = div_j_end

    hours_since_eruption    = int( time_since_eruption.total_seconds() / 3600 )
    minutes_since_eruption = int( ( time_since_eruption.total_seconds() % 3600 ) / 60 )
        
    title  = 'Estimated wavefront from Hunga Tonga\n'
    title += f'{hours_since_eruption:3d}:{minutes_since_eruption:02d} since the eruption'
    title += '[' + date_time.astimezone( timezone.utc ).strftime('%Y-%m-%d %H:%M:%S (UTC)') + ']\n'
    title += f'(Assumed travel speed = {travel_speed_m_s:.1f} m/s)'
    ax.set_title( title )

    if GENERATE_ANIMATION_GIF:
        print(f'Internally generated the map at {hours_since_eruption:3d}:{minutes_since_eruption:02d} since the eruption.')
        animation_images.append( fig2img( fig ) )
    else:
        fig_output_filename = OUTPUT_DIR / f'wavefront_simulation_{minutes_from_eruption:05d}.png'
        fig.savefig( fig_output_filename, bbox_inches='tight', pad_inches = 0 )
        print( f'Generated {fig_output_filename}' )
        
    fig.clear()


if GENERATE_ANIMATION_GIF:
    # Generate the cover page.
    image_size = animation_images[0].size
    cover_image = Image.new( 'RGB', image_size, 'black' )
    cover_draw  = ImageDraw.Draw( cover_image )
    cover_font  = ImageFont.truetype( findfont( FontProperties( family = ['monospace'] ) ),
                                      size = 24 )
    cover_text  = 'Estimated wavefront of the shockwave\nfrom Hunga Tonga'
    cover_draw.text( ( image_size[0] / 2, image_size[1] / 2 ),
                     cover_text,
                     fill = 'white',
                     font = cover_font,
                     anchor = 'mm',
                     align = 'center' )
    
    # Generate animation GIF.
    duration = [2000] + [100] * len( animation_images )
    gif_output_filename = OUTPUT_DIR / f'wavefront_simulation.gif'
    cover_image.save( gif_output_filename,
                      append_images = animation_images,
                      save_all = True,
                      duration = duration,
                      loop = 0 )
    print( f'Generated {gif_output_filename}' )
