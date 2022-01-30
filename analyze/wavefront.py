
#
# Python module that defines common functions and classes related to
# estimated wavefront from Hunga Tonga.
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

import numpy as np

from .common import *

class WavefrontLine:
    def __init__(self, travel_speed_m_s, color):
        self.travel_speed_m_s = travel_speed_m_s
        self.color            = color

def draw_wavefront(ax, distance, projection, wavefront_line):
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
                        color = wavefront_line.color )
        lines.append( line )
        div_j_start = div_j_end

    return lines
