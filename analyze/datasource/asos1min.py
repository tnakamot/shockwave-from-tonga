#
# Library to collect data from ASOS1MIN network through IEC (Iowa
# Environmental Mesonet).
#   https://mesonet.agron.iastate.edu/request/asos/1min.phtml
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

import csv
import urllib.request
from datetime import datetime, timedelta, timezone
from io import StringIO

from bs4 import BeautifulSoup
from tqdm import tqdm

from .base import BarometricPressureDataSourceBase
from ..database import BarometricPressureMonitoringStation, BarometricPressureRecord
from ..common import *

class BarometricPressureDataSourceAsosOneMinute( BarometricPressureDataSourceBase ):
    def __init__( self ):
        super( BarometricPressureDataSourceAsosOneMinute, self ).__init__( 'asos1min' )

    def _get_station_identifiers( self ):
        '''This method obtains a list of ASOS station identifiers in the ASOS1MIN network.

        Returns:
        --------
        list
            List of ASOS station identifiers (str)
        '''

        fp = urllib.request.urlopen( 'https://mesonet.agron.iastate.edu/plotting/asos/1station_1min.phtml' )
        soup = BeautifulSoup( fp, 'html.parser' )
        options = soup.select('select[name=station]')[0].find_all('option')
        return [ option['value'] for option in options ]

    def _get_station( self, station_identifier ):
        '''This method obtains the metadata of the specified ASOS station.

        Parameters:
        -----------
        station_identifier: str
            ASOS station identifier. E.g., CQT.

        Returns:
        --------
        AsosStation
            An instance of AsosStation where the station metadata is stored. E.g., station name.
        '''
        
        fp = urllib.request.urlopen( f'https://mesonet.agron.iastate.edu/sites/site.php?station={station_identifier}' )
        soup = BeautifulSoup( fp, 'html.parser' )
        tables = soup.find_all('table')

        name = ''
        network = ''
        latitude_deg = None
        longitude_deg = None
        elevation_m = None

        # TODO: make the following code more robust
        for table in tables:
            if table.tr.th.text != 'Station Identifier:':
                continue

            for row in table.find_all('tr'):
                if row.th.text == 'Station Name:':
                    name = row.td.text
                elif row.th.text == 'Network':
                    network = row.td.text
                elif row.th.text == 'Latitude:':
                    latitude_deg = float( row.td.text )
                elif row.th.text == 'Longitude:':
                    longitude_deg = float( row.td.text )
                elif row.th.text == 'Elevation [m]:':
                    elevation_m = float( row.td.text )
            break

        if latitude_deg is None or longitude_deg is None:
            return None

        return BarometricPressureMonitoringStation(
            data_source   = self.name,
            identifier    = station_identifier,
            latitude_deg  = latitude_deg,
            longitude_deg = longitude_deg,
            elevation_m   = elevation_m,
            name          = name,
        )

    def _get_barometric_pressure_records( self, station, start_date_time, end_date_time ):
        '''This method downloads 1-minute pressure records of the specified ASOS station in the specified period.

        Parameters:
        -----------
        station: AsosStation
            The ASOS station.

        start_date_time: datetime
            The beginning of the period for the data query.

        end_date: datetime
            The end of the period for the data query.

        Returns:
        --------
        list
             1-minute pressure records (AsosRecord)
        '''
        
        start_utc = start_date_time.astimezone( timezone.utc )
        end_utc   = end_date_time.astimezone( timezone.utc )

        query_params = {
            'station[]': [station.identifier] ,
            'tz'       : 'UTC',
            'year1'    : start_utc.year,
            'month1'   : start_utc.month,
            'day1'     : start_utc.day,
            'hour1'    : start_utc.hour,
            'minute1'  : start_utc.minute,
            'year2'    : end_utc.year,
            'month2'   : end_utc.month,
            'day2'     : end_utc.day,
            'hour2'    : end_utc.hour,
            'minute2'  : end_utc.minute,
            'vars[]'   : ['pres1'],
            'sample'   : '1min',
            'what'     : 'download',
            'delim'    : 'comma',
            'gis'      : 'no',
        }

        BASE_URL  = 'https://mesonet.agron.iastate.edu/request/asos/1min_dl.php'
        query_url = BASE_URL + '?' + urllib.parse.urlencode( query_params, doseq = True )

        fp = urllib.request.urlopen( query_url )
        csv_reader = csv.DictReader( StringIO ( fp.read().decode() ) )
        records = []
    
        for row in csv_reader:
            if row['station'] != station.identifier:
                continue

            if not row['valid(UTC)']:
                continue

            date_time = datetime.strptime( row['valid(UTC)'], '%Y-%m-%d %H:%M' )
            date_time = datetime( date_time.year, date_time.month,  date_time.day,
                                  date_time.hour, date_time.minute, tzinfo = timezone.utc )

            if not row['pres1']:
                continue
        
            pres1_inchHg = float(row['pres1'])
            pres1_hPa    = pres1_inchHg * 33.863886666667

            record = BarometricPressureRecord( station, date_time, pres1_hPa )
            records.append( record )

        return records
        

    def import_data( self, database ):
        START_TIME      = ERUPTION_TIME
        END_TIME        = ERUPTION_TIME + timedelta( hours = 216 )

        tqdm.write( 'Obtaining a list of ASOS stations in the ASOS1MIN network... ', end = '' )
        station_identifiers = self._get_station_identifiers()
        tqdm.write( f'{len(station_identifiers)} ASOS stations were found.\n' )

        total_records = 0

        for station_identifier in tqdm( station_identifiers ):
            tqdm.write( f'Downloading data of ASOS station {station_identifier}... ', end = '' )
            station = self._get_station( station_identifier )
            if station is None:
                continue
            station.insert_to( database )

            records = self._get_barometric_pressure_records( station, START_TIME, END_TIME )
            BarometricPressureRecord.insert_records( database, records )
            total_records += len( records )
            
            tqdm.write( f'imported { len( records ) } records into the database.' )
            
        tqdm.write( f'Imported { total_records } records into the database { str( database ) }.' )
