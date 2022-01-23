#!/usr/bin/env python3

#
# Program to collect barometric pressure data from ASOS1MIN network
# through IEC (Iowa Environmental Mesonet).
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
import urllib.parse
import sqlite3
from datetime import datetime, timedelta, timezone
from io import StringIO
from math import nan
from pathlib import Path

from bs4 import BeautifulSoup
from geopy.distance import geodesic
from tqdm import tqdm

from common import *

STATION_TABLE_NAME = 'station'
PRESSURE_DATA_TABLE_NAME = 'pressure_data'

class AsosStation:
    '''This class represents one ASOS station.'''

    def __init__(self, identifier, name, network, latitude_deg, longitude_deg, elevation_m):
        self.identifier    = identifier
        self.name          = name
        self.network       = network
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.elevation_m   = elevation_m
        self.distance_km   = geodesic( HUNGA_TONGA_COORD, ( latitude_deg, longitude_deg ) ).km

    def insert(self, sqlite3_connection):
        sqlite3_connection.execute( f'''
INSERT INTO {STATION_TABLE_NAME}
(id, name, network, latitude_deg, longitude_deg, elevation_m, distance_km)
VALUES
(?, ?, ?, ?, ?, ?, ?)
''',
        [self.identifier,
         self.name,
         self.network,
         self.latitude_deg,
         self.longitude_deg,
         self.elevation_m,
         self.distance_km] )
        sqlite3_connection.commit()
        
class AsosRecord:
    '''One record in the data from an ASOS station.'''

    def __init__(self, date_time, pressure_hPa):
        self.date_time    = date_time
        self.pressure_hPa = pressure_hPa
        
def get_asos1min_station_identifiers():
    '''This function obtains a list of ASOS station identifiers in the ASOS1MIN network.

    Returns:
    --------
    list
        List of ASOS station identifiers (str)
    '''

    fp = urllib.request.urlopen( 'https://mesonet.agron.iastate.edu/plotting/asos/1station_1min.phtml' )
    soup = BeautifulSoup( fp, 'html.parser' )
    options = soup.select('select[name=station]')[0].find_all('option')
    return [ option['value'] for option in options ]

def get_asos_station(station_identifier):
    '''This function obtains the metadata of the specified ASOS station.

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
    latitude_deg = nan
    longitude_deg = nan
    elevation_m = nan

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

    return AsosStation( station_identifier, name, network, latitude_deg, longitude_deg, elevation_m )

def get_asos1min_records(station, start_date_time, end_date_time):
    '''This function downloads 1-minute pressure records of the specified ASOS station in the specified period.

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
        
        records.append( AsosRecord( date_time, pres1_hPa ) )

    return records

def create_sqlite3_tables(sqlite3_connection):   
    print( f'Creating table "{STATION_TABLE_NAME}" ...' )
    sqlite3_connection.execute( f'DROP TABLE IF EXISTS {STATION_TABLE_NAME};' )
    create_station_table = f'''
CREATE TABLE {STATION_TABLE_NAME} (
    id TEXT PRIMARY KEY,
    name          TEXT,
    network       TEXT,
    latitude_deg  REAL,
    longitude_deg REAL,
    elevation_m   REAL,
    distance_km   REAL
);
'''
    sqlite3_connection.execute( create_station_table )

    print( f'Creating table "{PRESSURE_DATA_TABLE_NAME}" ...' )
    sqlite3_connection.execute( f'DROP TABLE IF EXISTS {PRESSURE_DATA_TABLE_NAME};' )

    create_pressure_data_table = f'''
CREATE TABLE {PRESSURE_DATA_TABLE_NAME} (
    station_id    TEXT,
    timestamp     TIMESTAMP,
    pressure_hPa  REAL,
    FOREIGN KEY(station_id) REFERENCES station(id),
    PRIMARY KEY (station_id, timestamp)
);
'''
    sqlite3_connection.execute( create_pressure_data_table )

def insert_records(sqlite3_connection, station, records):
    sqlite3_connection.executemany( f'''
INSERT INTO {PRESSURE_DATA_TABLE_NAME}
(station_id, timestamp, pressure_hPa)
VALUES
(?, ?, ?)
    ''', [ (station.identifier, record.date_time, record.pressure_hPa) for record in records ] )

    sqlite3_connection.commit()

def main():
    '''Main function'''

    START_TIME      = ERUPTION_TIME
    END_TIME        = ERUPTION_TIME + timedelta( hours = 168 )

    print( f'Opening SQLite3 database {ASOS1MIN_SQLITE3_DATABASE} ...' )
    sqlite3_connection = sqlite3.connect( ASOS1MIN_SQLITE3_DATABASE )

    create_sqlite3_tables( sqlite3_connection )    

    print( 'Obtaining a list of ASOS stations in the ASOS1MIN network... ', end='', flush=True )
    station_identifiers = get_asos1min_station_identifiers()
    print( f'{len(station_identifiers)} ASOS stations were found.\n' )

    total_records = 0
    
    for station_identifier in tqdm( station_identifiers ):
        tqdm.write( f'Downloading data of ASOS station {station_identifier}... ', end = '' )
        station = get_asos_station( station_identifier )
        station.insert( sqlite3_connection )
        
        records = get_asos1min_records( station, START_TIME, END_TIME )
        insert_records( sqlite3_connection, station, records )
        total_records += len( records )
        
        tqdm.write( f'imported { len( records ) } records into the SQLite3 database.' )

    tqdm.write( f'Imported { total_records } records into the SQLite3 database {ASOS1MIN_SQLITE3_DATABASE} (table = {PRESSURE_DATA_TABLE_NAME} )' )

    sqlite3_connection.close()
    
if __name__ == "__main__":
    main()
