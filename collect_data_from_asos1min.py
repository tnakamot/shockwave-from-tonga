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
from datetime import datetime, timedelta, timezone
from io import StringIO
from math import nan
from pathlib import Path

from bs4 import BeautifulSoup
from geopy.distance import geodesic
from influxdb import InfluxDBClient
from tqdm import tqdm

from common import *

INFLUX_DB_MEASUREMENT = 'asos1min' # This must match the one in plot_asos1min.py

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

def convert_asos_records_to_csv(station, records):
    '''This function outputs CSV string of the specified ASOS pressure records.

    Parameters:
    -----------
    station: AsosStation
        ASOS station where the records were obtained.
    
    records: list
        ASOS pressure data records.

    Returns:
    --------
    str
         CSV in string.
    '''

    s = f'''
Station Identifier       , {station.identifier}
Station Name             , {station.name}
Network                  , {station.network}
Latitude [deg]           , {station.latitude_deg}
Longitude [deg]          , {station.longitude_deg}
Elevation [m]            , {station.elevation_m}

Date & Time              , Prs [inch]
'''
    
    for record in records:
        s += f'{record.date_time.isoformat()}, {record.sensor1_inch:11.6f}\n'
        
    return s

def record_to_influxdb_point(station, record):
    return {
        'measurement': INFLUX_DB_MEASUREMENT,
        'tags': {
            'station_id'   : station.identifier,
            'latitude_deg' : station.latitude_deg,
            'longitude_deg': station.longitude_deg,
            'distance_km'  : station.distance_km,
        },
        "time": record.date_time,
        "fields": {
            'pressure_hPa': record.pressure_hPa
        }
    }
    pass

def main():
    '''Main function'''

    START_TIME      = ERUPTION_TIME
    END_TIME        = ERUPTION_TIME + timedelta( hours = 144 )
    DATA_OUTPUT_DIR = Path( 'data_asos1min' )

    DATA_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

    print( f'Connecting InfluxDB on {INFLUX_DB_HOST}:{INFLUX_DB_PORT} ...' )
    influx_client = InfluxDBClient( INFLUX_DB_HOST , INFLUX_DB_PORT )
    
    print( f'Creating InfluxDB database "{INFLUX_DB_NAME}" ...' )
    influx_client.create_database( INFLUX_DB_NAME )
    influx_client.switch_database( INFLUX_DB_NAME )

    print( 'Obtaining a list of ASOS stations in the ASOS1MIN network... ', end='', flush=True )
    station_identifiers = get_asos1min_station_identifiers()
    print( f'{len(station_identifiers)} ASOS stations were found.\n' )

    total_records = 0
    
    for station_identifier in tqdm( station_identifiers ):
        tqdm.write( f'Downloading data of ASOS station {station_identifier}... ', end = '' )
        station = get_asos_station( station_identifier )
        records = get_asos1min_records( station, START_TIME, END_TIME )
        total_records += len( records )
        
        influx_client.write_points( [ record_to_influxdb_point( station, record)
                                      for record in records ] )

        tqdm.write( f'imported { len( records ) } data points into InfluxDB (measurement = { INFLUX_DB_MEASUREMENT })' )

        # Do not use CSV file anymore.
        #
        # csv_str = convert_asos_records_to_csv(station, records)
        # output_file = DATA_OUTPUT_DIR / f'{station.identifier}.csv'
        # f = open( output_file, 'w', encoding = 'utf-8' )
        # f.write( csv_str )
        # f.close()

        # tqdm.write( f'saved in {output_file}' )

    tqdm.write( f'Imported { total_records } data points in total into InfluxDB (measurement = { INFLUX_DB_MEASUREMENT })' )
    
if __name__ == "__main__":
    main()
