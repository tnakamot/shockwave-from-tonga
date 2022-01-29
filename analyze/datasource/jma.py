#
# Library to collect data from Japan Meteorological Agency.
#   https://www.jma.go.jp/
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

import re
import urllib.request
from datetime import datetime, timedelta, timezone

from bs4 import BeautifulSoup
from tqdm import tqdm

from .base import BarometricPressureDataSourceBase
from ..database import BarometricPressureMonitoringStation, BarometricPressureRecord
from ..common import *

class JmaArea:
    '''This class represents one area of JMA.'''
    
    def __init__( self, name, prec_no ):
        self.name    = name    # Area name   (str)
        self.prec_no = prec_no # Area number (int)

class JmaStation( BarometricPressureMonitoringStation ):
    '''This class represents one weather station of JMA.'''
    
    def __init__( self, data_source, area, block_no, name, latitude_deg, longitude_deg, elevation_m ):
        super( JmaStation, self ).__init__(
            data_source   = data_source,
            identifier    = str(block_no),
            latitude_deg  = latitude_deg,
            longitude_deg = longitude_deg,
            elevation_m   = elevation_m,
            name          = name,
        )
        self.area          = area     # Instance of Jma Area
        self.block_no      = block_no # Block number (int)

class BarometricPressureDataSourceJMA( BarometricPressureDataSourceBase ):
    def __init__( self ):
        super( BarometricPressureDataSourceJMA, self ).__init__( 'jma' )

    def _get_areas( self ):
        '''This function obtains JMA area information from the JMA's web page.

        Returns:
        --------
        list
            List of JmaArea instances.
        '''
        fp = urllib.request.urlopen( 'https://www.data.jma.go.jp/obd/stats/etrn/select/prefecture00.php' )
        soup = BeautifulSoup( fp, 'html.parser' )
        return [ JmaArea( name    = area_tag['alt'],
                          prec_no = int(re.search(r'prec_no=([0-9]+)', area_tag['href']).group(1)) )
                 for area_tag in soup.find_all('area') ]

    def _get_stations( self, area ):
        '''This function obtains information of all JMA weather stations in the specified JMA area.

        Parameters:
        -----------
        area: JmaArea
            The JMA area in which you want to search for JMA stations.

        Returns:
        --------
        list
            List of JmaStation instances.
        '''
    
        fp = urllib.request.urlopen( f'https://www.data.jma.go.jp/obd/stats/etrn/select/prefecture.php?prec_no={area.prec_no:02d}' )
        soup = BeautifulSoup(fp, 'html.parser')
        stations = []
        block_nos = {}

        for area_tag in soup.find_all('area'):
            if 'onmouseover' not in area_tag.attrs:
                continue
        
            args = re.search(r'viewPoint\((.*)\)', area_tag['onmouseover']).group(1).split(',')
            args = [arg.strip("'") for arg in args]

            if args[0] != 's':
                continue
            
            block_no = int( args[1] )
            if block_no in block_nos:
                continue
                
            station_name  = args[2]
            latitude_deg  = float( args[4] ) + ( float( args[5] ) / 60 )
            longitude_deg = float( args[6] ) + ( float( args[7] ) / 60 )
            elevation_m   = float( args[8] )
            
            station = JmaStation(
                data_source   = self.name,
                area          = area,
                block_no      = block_no,
                name          = station_name,
                latitude_deg  = latitude_deg,
                longitude_deg = longitude_deg,
                elevation_m   = elevation_m,
            )
            stations.append( station )

            block_nos[block_no] = True
        return stations

    def _get_ten_minutes_data_per_day( self, station, date_jst ):
        '''This function downloads 10-minute barometric pressure data of the specified weather station on the specified day.

        Parameters:
        -----------
        station: JmaStation
            The JMA weather station.

        date_jst: datetime
            The day on which you want to obtain 10-minute pressure data.

        Returns:
        --------
        list
            List of 10-minute barometric pressure records.
        '''
    
        prec_no  = station.area.prec_no
        block_no = station.block_no
        year     = date_jst.year
        month    = date_jst.month
        day      = date_jst.day

        records = []

        url = f'https://www.data.jma.go.jp/obd/stats/etrn/view/10min_s1.php?prec_no={prec_no:02d}&block_no={block_no:04d}&year={year:04d}&month={month:02d}&day={day:02d}&view='
        fp = urllib.request.urlopen( url )
        soup = BeautifulSoup( fp, 'html.parser' )

        data_table = soup.find( id = 'tablefix1' )
        data_table_rows = data_table.find_all( 'tr' )[2:]
        for row in data_table_rows:
            cells = row.find_all( 'td' )
            hour, minute = [ int( s ) for s in cells[0].string.split(':') ]
        
            if hour == 24:
                date_time = datetime( year, month, day, 0, 0, tzinfo = TZ_JST ) + timedelta( days = 1 )
            else:
                date_time = datetime( year, month, day, hour, minute, tzinfo = TZ_JST )

            try:
                pressure_hPa = float( cells[1].string )
            except:
                continue

            record = BarometricPressureRecord( station, date_time, pressure_hPa )
            records.append( record )
            
        return records

    def _get_ten_minutes_data( self, station, start_date, end_date ):
        '''This function downloads 10-minute pressure data of the specified weather station on the specified days.

        Parameters:
        -----------
        station: JmaStation
            The JMA weather station.

        start_date: datetime
            The first day for the data query.

        end_date: datetime
            The last day for the data query.

        Returns:
        --------
        list
            List of 10-minute barometric pressure records.
        '''

        current_date_jst = start_date.astimezone( TZ_JST )
        records = []
    
        while current_date_jst <= end_date:
            records.extend( self._get_ten_minutes_data_per_day( station, current_date_jst ) )
            current_date_jst += timedelta( days = 1 )
                            
        return records
    
    def import_data( self, database ):
        ERUPTION_TIME_JST = ERUPTION_TIME.astimezone( TZ_JST )
        START_DATE = datetime( ERUPTION_TIME_JST.year,
                               ERUPTION_TIME_JST.month,
                               ERUPTION_TIME_JST.day,
                               tzinfo = TZ_JST )
        END_DATE   = datetime( ERUPTION_TIME_JST.year,
                               ERUPTION_TIME_JST.month,
                               ERUPTION_TIME_JST.day + 4,
                               tzinfo = TZ_JST )

        tqdm.write( 'Obtaining JMA area information... ', end = '' )
        areas = self._get_areas()
        tqdm.write( f'{len(areas)} JMA areas were found.\n' )

        stations = []
        for area in tqdm( areas ):
            tqdm.write( f'Obtaining JMA weather stations in area {area.name}...' )
            stations.extend( self._get_stations( area ) )
        tqdm.write( f'\n{len(stations)} JMA weather stations were found.\n' )
        
        total_records = 0
        for station in tqdm( stations ):
            station_name_with_padding = station.name + ' ' * ( (8 - len(station.name)) * 2 )
            tqdm.write( f'Downloading data from weather station {station_name_with_padding}... ', end = '' )
            station.insert_to( database )

            records = self._get_ten_minutes_data( station, START_DATE, END_DATE )
            BarometricPressureRecord.insert_records( database, records )
            tqdm.write( f'imported { len( records ) } records into the database.' )
            
            total_records += len( records )
            
        tqdm.write( f'Imported { total_records } records into the database { str( database ) }.' )
