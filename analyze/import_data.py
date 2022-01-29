#!/usr/bin/env python3

#
# Program to import data from various data sources.
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

from .datasource.asos1min import BarometricPressureDataSourceAsosOneMinute
from .datasource.jma import BarometricPressureDataSourceJMA
from .database import AnalysisDatabase

DATA_SOURCES = {
    'asos1min': BarometricPressureDataSourceAsosOneMinute(),
    'jma'     : BarometricPressureDataSourceJMA(),
}

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description = 'This program imports data related to shockwave from Hunga Tonga from various data sources and store them in the database.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        '--db',
        default = AnalysisDatabase.DEFAULT_DBFILE_PATH,
        help = 'Path to the SQLite3 database file.',
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
    selected_data_source = DATA_SOURCES[ args.source ]
    selected_data_source.import_data( database )

    database.close()

if __name__ == "__main__":
    main()
