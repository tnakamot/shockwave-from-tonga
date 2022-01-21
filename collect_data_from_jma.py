#!/usr/bin/env python3

#
# Program to collect barometric pressure data from Japan Meteorological Agency.
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
from math import nan
from pathlib import Path

from bs4 import BeautifulSoup
from tqdm import tqdm

TZ_JST = timezone( timedelta( hours = +9 ), name = 'JST' ) # Japan Standard Time

class JmaArea:
    '''This class represents one area of JMA.'''
    
    def __init__(self, name, prec_no):
        self.name    = name    # Area name   (str)
        self.prec_no = prec_no # Area number (int)

class JmaStation:
    '''This class represents one weather station of JMA.'''
    
    def __init__(self, area, block_no, name, latitude_deg, longitude_deg, height_m):
        self.area          = area            # Instance of Jma Area
        self.block_no      = block_no        # Block number (int)
        self.name          = name            # Name of the weather station (str)
        self.latitude_deg  = latitude_deg    # Latitude in degree (float)
        self.longitude_deg = longitude_deg   # Longitude in degree (float)
        self.height_m      = height_m        # Height in meter (float)

class JmaTenMinutesData:
    '''An instance of this class contains 10-minute pressure data of one JMA weather station.'''
    
    def __init__(self,
                 station,
                 date_time,
                 pressure_hPa,
                 sea_level_pressure_hPa):
        self.station                = station
        self.date_time              = date_time
        self.pressure_hPa           = pressure_hPa
        self.sea_level_pressure_hPa = sea_level_pressure_hPa
        
def get_jma_areas():
    '''This function obtains JMA area information from the JMA's web page.

    Returns:
    --------
    list
        List of JmaArea instances.
    '''
    fp = urllib.request.urlopen('https://www.data.jma.go.jp/obd/stats/etrn/select/prefecture00.php')
    soup = BeautifulSoup(fp, 'html.parser')
    return [ JmaArea( name    = area_tag['alt'],
                      prec_no = int(re.search(r'prec_no=([0-9]+)', area_tag['href']).group(1)) )
             for area_tag in soup.find_all('area') ]

def get_jma_stations(area):
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
    
    fp = urllib.request.urlopen(f'https://www.data.jma.go.jp/obd/stats/etrn/select/prefecture.php?prec_no={area.prec_no:02d}')
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
        height_m      = float( args[8] )
            
        station = JmaStation( area          = area,
                              block_no      = block_no,
                              name          = station_name,
                              latitude_deg  = latitude_deg,
                              longitude_deg = longitude_deg,
                              height_m      = height_m )
        stations.append( station )

        block_nos[block_no] = True
        
    return stations

def get_jma_ten_minutes_data_per_day(station, date_jst):
    '''This function downloads 10-minute pressure data of the specified weather station on the specified day.

    Parameters:
    -----------
    station: JmaStation
        The JMA weather station.

    date_jst: datetime
        The day on which you want to obtain 10-minute pressure data.

    Returns:
    --------
    date_time: list
        List of the timestamp (datetime).

    pressure_hPa: list
        List of observed pressure in hPa (float).

    sea_level_pressure_hPa: list
        List of observed sea-level equivalent pressure in hPa (float).
    '''
    
    prec_no  = station.area.prec_no
    block_no = station.block_no
    year     = date_jst.year
    month    = date_jst.month
    day      = date_jst.day

    date_time    = []
    pressure_hPa = []
    sea_level_pressure_hPa = []

    url = f'https://www.data.jma.go.jp/obd/stats/etrn/view/10min_s1.php?prec_no={prec_no:02d}&block_no={block_no:04d}&year={year:04d}&month={month:02d}&day={day:02d}&view='
    fp = urllib.request.urlopen( url )
    soup = BeautifulSoup( fp, 'html.parser' )

    data_table = soup.find( id = 'tablefix1' )
    data_table_rows = data_table.find_all( 'tr' )[2:]
    for row in data_table_rows:
        cells = row.find_all( 'td' )
        hour, minute = [ int( s ) for s in cells[0].string.split(':') ]
        
        if hour == 24:
            date_time.append( datetime( year, month, day, 0, 0, tzinfo = TZ_JST ) + timedelta( days = 1 ) )
        else:
            date_time.append( datetime( year, month, day, hour, minute, tzinfo = TZ_JST ) )

        try:
            p = float( cells[1].string )
        except:
            p = nan
        pressure_hPa.append( p )

        try:
            p = float( cells[2].string )
        except:
            p = nan
        sea_level_pressure_hPa.append( p )

    return date_time, pressure_hPa, sea_level_pressure_hPa

def get_jma_ten_minutes_data(station, start_date, end_date):
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
    JmaTenMinutesData
         10-minute pressure data of the specified weather station on the specified days.
    '''

    current_date_jst = start_date.astimezone( TZ_JST )
    date_time = []
    pressure_hPa = []
    sea_level_pressure_hPa = []
    
    while current_date_jst <= end_date:
        r1, r2, r3 = get_jma_ten_minutes_data_per_day( station, current_date_jst )
        date_time.extend( r1 )
        pressure_hPa.extend( r2 )
        sea_level_pressure_hPa.extend( r3 )
    
        current_date_jst += timedelta( days = 1 ) # Go to the next day

    return JmaTenMinutesData( station, date_time, pressure_hPa, sea_level_pressure_hPa )

def convert_ten_minutes_data_to_csv(data):
    '''This function outputs CSV string of the specified 10-minute pressure data.

    Parameters:
    -----------
    data: JmaTenMinuteData
        10-minute pressure data you want to save in the CSV format.

    Returns:
    --------
    str
         CSV in string.
    '''
    
    station = data.station
    s = f'''
prec_no                  , {station.area.prec_no:02d}
block_no                 , {station.block_no:02d}
Station name             , {station.name}
Latitude [deg]           , {station.latitude_deg:8.3f}
Longitude [deg]          , {station.longitude_deg:8.3f}
Height [m]               , {station.height_m:6.1f}

Date & Time              , Prs [hPa], Sea Level Prs [hPa]
'''
    for i, t in enumerate( data.date_time ):
        s += f'{t.isoformat()}, {data.pressure_hPa[i]:9.1f}, {data.sea_level_pressure_hPa[i]:9.1f}\n'
    return s


def main():
    '''Main function'''

    START_DATE      = datetime( 2022, 1, 15, tzinfo = TZ_JST ) # Analysis start date in JST (inclusive)
    END_DATE        = datetime( 2022, 1, 19, tzinfo = TZ_JST ) # Analysis end date in JST (inclusive)
    DATA_OUTPUT_DIR = Path( 'data_jma' )

    DATA_OUTPUT_DIR.mkdir( parents = True, exist_ok = True )

    print( 'Obtaining JMA area information... ', end='', flush=True )
    areas = get_jma_areas()
    print( f'{len(areas)} JMA areas were found.\n' )

    stations = []
    for area in tqdm( areas ):
        tqdm.write( f'Obtaining JMA weather stations in area {area.name}...' )
        stations.extend( get_jma_stations( area ) )

    print( f'\n{len(stations)} JMA weather stations were found.\n' )

    for station in tqdm( stations ):
        station_name_with_padding = station.name + ' ' * ( (8 - len(station.name)) * 2 )
        tqdm.write( f'Downloading data from weather station {station_name_with_padding}... ', end = '' )
        
        data = get_jma_ten_minutes_data( station, START_DATE, END_DATE )
        output_file = DATA_OUTPUT_DIR / f'{station.area.prec_no:02d}_{station.block_no:04d}.csv'
        f = open( output_file, 'w', encoding = 'utf-8' )
        f.write( convert_ten_minutes_data_to_csv( data ) )
        f.close()

        tqdm.write( f'saved in {output_file}' )

if __name__ == "__main__":
    main()
