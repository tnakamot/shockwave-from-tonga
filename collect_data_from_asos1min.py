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
from tqdm import tqdm

from common import *

class AsosStation:
    '''This class represents one ASOS station.'''

    def __init__(self, identifier, name, network, latitude_deg, longitude_deg, elevation_m):
        self.identifier    = identifier
        self.name          = name
        self.network       = network
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.elevation_m   = elevation_m

class AsosRecord:
    '''One record in the data from an ASOS station.'''

    def __init__(self, date_time, sensor1_inch, sensor2_inch, sensor3_inch):
        self.date_time    = date_time
        self.sensor1_inch = sensor1_inch
        self.sensor2_inch = sensor2_inch
        self.sensor3_inch = sensor3_inch
        
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
        'vars[]'   : ['pres1', 'pres2', 'pres3'],
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

        pres1 = float(row['pres1']) if row['pres1'] else nan
        pres2 = float(row['pres2']) if row['pres2'] else nan
        pres3 = float(row['pres3']) if row['pres3'] else nan

        records.append( AsosRecord( date_time, pres1, pres2, pres3 ) )

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

Date & Time              , Prs1 [inch], Prs2 [inch], Prs3 [inch]
'''
    
    for record in records:
        s += f'{record.date_time.isoformat()}, {record.sensor1_inch:11.6f}, {record.sensor2_inch:11.6f}, {record.sensor3_inch:11.6f}\n'
        
    return s

def main():
    '''Main function'''

    START_TIME      = ERUPTION_TIME
    END_TIME        = ERUPTION_TIME + timedelta( hours = 144 )
    DATA_OUTPUT_DIR = Path( 'data_asos1min' )

    DATA_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

    print( 'Obtaining a list of ASOS stations in the ASOS1MIN network... ', end='', flush=True )
    station_identifiers = get_asos1min_station_identifiers()
    print( f'{len(station_identifiers)} ASOS stations were found.\n' )
    
    for station_identifier in tqdm( station_identifiers ):
        tqdm.write( f'Downloading data of ASOS station {station_identifier}... ', end = '' )
        station = get_asos_station( station_identifier )
        records = get_asos1min_records( station, START_TIME, END_TIME )
        csv_str = convert_asos_records_to_csv(station, records)

        output_file = DATA_OUTPUT_DIR / f'{station.identifier}.csv'
        f = open( output_file, 'w', encoding = 'utf-8' )
        f.write( csv_str )
        f.close()

        tqdm.write( f'saved in {output_file}' )
    
if __name__ == "__main__":
    main()
