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
from .database import AnalysisDatabase, BarometricPressureMonitoringStation
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

        return ( ERUPTION_TIME + min( time_since_eruption ),
                 ERUPTION_TIME + max( time_since_eruption ) )

DATA_SOURCES = {
    'asos1min': ( BarometricPressureDataSourceAsosOneMinute(),
                  MapRange( Point( latitude = 25.0, longitude = - 65.0 ),
                            Point( latitude = 50.0, longitude = -125.0 ) ) ),
    'jma'     : ( BarometricPressureDataSourceJMA(),
                  MapRange( Point( latitude = 24.0, longitude = 120.0 ),
                            Point( latitude = 46.0, longitude = 160.0 ) ) ),
}

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

    animation_data = AnimationData( dump_mode = (dump_dir is not None) )

    # TODO: adjust accordingly
    MAX_PRESSURE_DIFF_HPA_MINUTE = 0.5

    projection = ccrs.PlateCarree()

    stations = BarometricPressureMonitoringStation.get_all_stations_in_datasource( database, data_source )
    tqdm.write( f'{len(stations)} stations were found' )

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
