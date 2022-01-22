
#
# Python module that defines common constants, functions and classes that are
# useful for analysing and visualizing the shockwave from Hunga Tonga.
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

import contextlib
import io
from datetime import datetime, timedelta, timezone

import cartopy.feature as cfea
import cv2
import numpy as np
from geopy.distance import geodesic
from geopy.point import Point
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties, findfont
from PIL import Image, ImageDraw, ImageFont

# Japan Standard Time
TZ_JST = timezone( timedelta( hours = +9 ), name = 'JST' ) # Japan Standard Time

OCEAN_COLOR = '#00FFFF'
LAND_COLOR  = '#32CD32'

# Coordinate of Hunga Tonga.
HUNGA_TONGA_COORD   = Point( latitude = -20.536 , longitude = -175.382 )

# Coordinate of antipode of Hunga Tonga.
ANTIPODE_COORD      = Point( latitude  = - HUNGA_TONGA_COORD.latitude,
                             longitude = 180 + HUNGA_TONGA_COORD.longitude )

# Earth circumference.
EARTH_CIRCUMFERENCE = geodesic( HUNGA_TONGA_COORD, ANTIPODE_COORD ) * 2

# Hunga Tonga eruption time estimated from the satellite image of Himawari 8
#  https://himawari.asia/himawari8-image.htm?sI=D531106&sClC=ffff00&sTA=true&sTAT=TY&sS=6&sNx=3&sNy=2&sL=-169.171875&sT=-426.8125&wW=1920&wH=969&sD=1642219800000
ERUPTION_TIME       = datetime( 2022, 1, 15, 13, 10, tzinfo = TZ_JST )

def ordinal(i):
    if 11 <= (i % 100) <= 13:
        return str(i) + 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd'] + ['th'] * 6
        return str(i) + suffix[i % 10]

def fig2img(fig):
    '''Convert the given matplotlib Figure into Pillow Image.

    Parameters:
    -----------
    fig: matplotlib.figure.Figure
        Figure to convert.

    Returns:
    --------
    PIL.Image
        Conversion result.
    '''
    buf = io.BytesIO()
    fig.savefig( buf, bbox_inches = 'tight', pad_inches = 0 )
    buf.seek( 0 )
    return Image.open( buf )

def add_world_map(fig, projection):
    ax = fig.add_subplot( 1, 1, 1, projection = projection )
    ax.set_extent( ( -180, 180, -90, 90 ), projection )
    ax.add_feature( cfea.OCEAN, color = OCEAN_COLOR )
    ax.add_feature( cfea.LAND,  color = LAND_COLOR )
    return ax

def draw_wavefront(ax, distance, projection, color = 'black'):
    bearings = np.linspace(-180, 180, 360)
    wavefront_points = [ distance.destination( point = HUNGA_TONGA_COORD, bearing = b ) for b in bearings ]

    wavefront_latitude_deg  = [ p[0] for p in wavefront_points ]
    wavefront_longitude_deg = [ p[1] for p in wavefront_points ]

    div_js = np.where( np.abs( np.diff( wavefront_longitude_deg ) ) > 180 )[0] + 1
    div_js = np.append( div_js, len( wavefront_longitude_deg ) )
    div_j_start = 0
    lines = []
    for div_j_end in div_js:
        line = ax.plot( wavefront_longitude_deg[div_j_start:div_j_end],
                        wavefront_latitude_deg[div_j_start:div_j_end],
                        transform = projection,
                        color = color )
        lines.append( line )
        div_j_start = div_j_end

    return lines


class AnimationFrame:
    '''Represents one image in the animation data.'''
    def __init__(self, image, duration_ms):
        self.image       = image
        self.duration_ms = duration_ms

class AnimationData:
    '''Represents one animation data.'''

    def __init__(self, dump_mode = True):
        self.frames = []
        self.cover_frame = None
        self.dump_mode   = dump_mode
        
    def append(self, image, duration_ms, dump_path = None):
        if isinstance(image, Figure):
            image = fig2img( image )
            
        self.frames.append( AnimationFrame( image, duration_ms ) )
        if self.dump_mode and dump_path:
            image.save( dump_path )

    def add_cover(self, text, fontsize = 24, duration_ms = 2000):
        image_size  = self.frames[0].image.size
        cover_image = Image.new( 'RGB', image_size, 'black' )
        cover_draw  = ImageDraw.Draw( cover_image )
        cover_font  = ImageFont.truetype( findfont( FontProperties( family = ['monospace'] ) ),
                                          size = fontsize )
        cover_draw.text( ( image_size[0] / 2, image_size[1] / 2 ),
                         text,
                         fill   = 'white',
                         font   = cover_font,
                         anchor = 'mm',
                         align  = 'center' )
        
        self.cover_frame = AnimationFrame( cover_image, duration_ms )
        
    def save_gif(self, path, loop = True):
        frames = [self.cover_frame] if self.cover_frame else []
        frames = frames + self.frames

        frames[0].image.save( path,
                              append_images = [ frame.image for frame in frames[1:] ],
                              save_all = True,
                              duration = [ frame.duration_ms for frame in frames ],
                              loop = 0 if loop else 1 )
    def save_mp4(self, path):
        frames = [self.cover_frame] if self.cover_frame else []
        frames = frames + self.frames

        one_frame_duration_ms = min( [ frame.duration_ms for frame in frames ] )
        fps = 1000.0 / one_frame_duration_ms

        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video = cv2.VideoWriter( str( path ), fourcc, fps, self.frames[0].image.size )

        for frame in frames:
            cv_image = cv2.cvtColor( np.array( frame.image ), cv2.COLOR_RGB2BGR )
            
            for i in range( int( frame.duration_ms / one_frame_duration_ms ) ):
                video.write( cv_image )
            
        video.release()
