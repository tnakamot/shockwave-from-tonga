#
# Base classes to create and access barometric pressure database using SQLite3.
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

import os
import sqlite3
from pathlib import Path

from geopy.distance import geodesic
from tqdm import tqdm

from .common import *

class AnalysisDatabase:
    '''Represents the database for analysis.'''

    DEFAULT_DBFILE_PATH = Path( os.path.dirname( os.path.realpath( __file__ ) ) ).parent / 'database.sqlite3'

    def __init__( self, dbfile_path = DEFAULT_DBFILE_PATH ):
        '''
        Parameters:
        -----------

        dbfile_path: str or Path
            The file path of SQLite3 database file.
        '''
        self.dbfile_path = dbfile_path

        tqdm.write( f'Opening SQLite3 database {self.dbfile_path} ...' )
        self._connection = sqlite3.connect( self.dbfile_path )

        BarometricPressureMonitoringStation.create_table( self )
        BarometricPressureRecord.create_table( self )

    def execute( self, sql_statement, parameters = None ):
        if parameters is None:
            self._connection.execute( sql_statement )
        else:
            self._connection.execute( sql_statement, parameters )
        self._connection.commit()

    def executemany( self, sql_statement, parameters ):
        self._connection.executemany( sql_statement, parameters )
        self._connection.commit()

    def cursor( self ):
        return self._connection.cursor()

    def table_exists( self, table_name ):
        c = self._connection.cursor()
        c.execute( f'''
SELECT
    name
FROM
    sqlite_master
WHERE
    type = 'table' AND
    name = '{table_name}'
''' )
        return c.fetchone() is not None

    def close( self ):
        self._connection.close()

    def __str__( self ):
        return str( self.dbfile_path )
       
class BarometricPressureMonitoringStation:
    '''This class represents one barometric pressure monitoring station.'''

    TABLE_NAME = 'barometric_pressure_monitor_station'

    def __init__(
            self,
            data_source,
            identifier,
            latitude_deg,
            longitude_deg,
            elevation_m = None,
            name = None,
    ):
        '''
        Parameters:
        -----------
        data_source: str
            Unique name of the data source (e.g., JMA, ASOS1MIN).

        identifier: str
            Identifier of the monitoring station. It must be unique within the same data source.

        latitude_deg: float
            Latitude of the monitoring station in degrees.

        longitude_deg: float
            Longitude of the monitoring station in degrees.

        elevation_m: float, optional
            Elevation of the monitoring station in meters.

        name: str, optional
            Name of the monitoring station.
        '''

        self.data_source   = data_source
        self.identifier    = identifier
        self.name          = name
        self.latitude_deg  = latitude_deg
        self.longitude_deg = longitude_deg
        self.elevation_m   = elevation_m
        self.distance_km   = geodesic( HUNGA_TONGA_COORD, ( latitude_deg, longitude_deg ) ).km

    @classmethod
    def create_table( cls, database ):
        '''Create the table of monitoring stations in the given SQLite3 database if not exists.'''

        if database.table_exists( cls.TABLE_NAME ):
            return
        
        tqdm.write( f'Creating table "{cls.TABLE_NAME}" ...' )
        
        sql_statement = f'''
CREATE TABLE IF NOT EXISTS {cls.TABLE_NAME} (
    data_source   TEXT NOT NULL,
    identifier    TEXT NOT NULL,
    latitude_deg  REAL NOT NULL,
    longitude_deg REAL NOT NULL,
    distance_km   REAL NOT NULL,
    elevation_m   REAL,
    name          TEXT,

    PRIMARY KEY (data_source, identifier)
);
'''
        database.execute( sql_statement )

    def insert_to( self, database ):
        database.execute( f'''
INSERT OR IGNORE INTO {self.TABLE_NAME}
(data_source, identifier, latitude_deg, longitude_deg, distance_km, elevation_m, name)
VALUES
(?, ?, ?, ?, ?, ?, ?)
''',
        [
            self.data_source,
            self.identifier,
            self.latitude_deg,
            self.longitude_deg,
            self.distance_km,
            self.elevation_m,
            self.name
        ] )

    @classmethod
    def _get_all( cls, database, data_source = None ):
        cursor = database.cursor()

        sql_statement = f'''
SELECT 
    data_source,
    identifier,
    latitude_deg,
    longitude_deg,
    distance_km,
    name
FROM
    {cls.TABLE_NAME}
'''

        # TODO: avoid SQL injection
        if data_source:
            sql_statement += f'''
WHERE
    data_source = "{data_source.name}"
'''
            
        sql_statement += '''
ORDER BY
    distance_km
'''

        cursor.execute( sql_statement )
            
        rows = cursor.fetchall()
        return [
            BarometricPressureMonitoringStation(
                data_source   = fields[0],
                identifier    = fields[1],
                latitude_deg  = fields[2],
                longitude_deg = fields[3],
                elevation_m   = fields[4],
                name          = fields[5],
            ) for fields in rows ]

    @classmethod
    def get_all_stations( cls, database ):
        return cls._get_all( database )

    @classmethod
    def get_all_stations_in_datasource( cls, database, data_source ):
        return cls._get_all( database, data_source )

class BarometricPressureMonitoringStationData:
    '''This class handles database access to the pressure data of the specified station.'''

    def __init__( self, database, station, start_time, end_time ):
        self._database   = database
        self._station    = station
        self._start_time = start_time.astimezone( timezone.utc )
        self._end_time   =   end_time.astimezone( timezone.utc )
        self._load()
        
    def _load( self ):
        cursor = self._database.cursor()
        cursor.execute( f'''
SELECT
    timestamp,
    pressure_hPa
FROM
    {BarometricPressureRecord.TABLE_NAME}
WHERE
    data_source = ? AND
    station_id  = ? AND
    timestamp  >= ? AND
    timestamp  <  ?
ORDER BY
    timestamp ASC
        ''', ( self._station.data_source,
               self._station.identifier,
               self._start_time,
               self._end_time ) )
        rows = cursor.fetchall()
        self._timestamps   = [ datetime.fromisoformat( row[0] ) for row in rows ]
        self._pressure_hPa = [ row[1] for row in rows ]

    def hours_since_eruption_range( self ):
        return ( ( self._start_time - ERUPTION_TIME ).total_seconds() / 3600,
                 ( self._end_time   - ERUPTION_TIME ).total_seconds() / 3600 )
    
    def is_empty( self ):
        '''Returns True if there is no pressure data or only one record in the specified date and time range.'''
        return len( self._timestamps ) <= 1

    def _get_pressure_hPa( self, interpolate = None ):
        minutes_since_eruption = np.array( [ ( t - ERUPTION_TIME ).total_seconds() / 60
                                             for t in self._timestamps ] )
        if not interpolate:
            return self._pressure_hPa, minutes_since_eruption

        start_minutes = int( ( self._start_time - ERUPTION_TIME ).total_seconds() / 60 ) + 1
        end_minutes   = int( ( self._end_time   - ERUPTION_TIME ).total_seconds() / 60 )
        interpolated_minutes_since_eruption = np.arange( start_minutes, end_minutes, 1 )
        f_interpolate_pressure_hPa = interp1d( minutes_since_eruption,
                                               self._pressure_hPa,
                                               kind = interpolate,
                                               copy = False,
                                               assume_sorted = True )
        extrapolate_start_n = int(minutes_since_eruption[0] - start_minutes)
        extrapolate_end_n   = int(end_minutes - minutes_since_eruption[-1]) 
        interp_start_i = extrapolate_start_n
        interp_end_i   = len( interpolated_minutes_since_eruption ) - extrapolate_end_n
        
        interpolated_pressure_hPa = np.concatenate( (
            np.full( extrapolate_start_n, self._pressure_hPa[0] ),
            f_interpolate_pressure_hPa( interpolated_minutes_since_eruption[ interp_start_i : interp_end_i ] ),
            np.full( extrapolate_end_n, self._pressure_hPa[-1] )
        ) )
        
        return interpolated_pressure_hPa, interpolated_minutes_since_eruption

    def _pressure_diff_hPa_minute_interp( self, interpolate ):
        pressure_hPa, minutes_since_eruption = self._get_pressure_hPa( interpolate )
        _pressure_diff_hPa_minute = np.diff( pressure_hPa )
        return _pressure_diff_hPa_minute, minutes_since_eruption[1:]
    
    def pressure_diff_hPa_minute( self, interpolate = None ):
        if interpolate:
            return self._pressure_diff_hPa_minute_interp( interpolate )

        previous_timestamp = None
        previous_pressure_hPa = None
        minutes_since_eruption = []
        _pressure_diff_hPa_minute = []
        for i, current_timestamp in enumerate( self._timestamps ):
            current_pressure_hPa = self._pressure_hPa[ i ]

            if previous_timestamp:
                seconds_since_eruption = ( current_timestamp - ERUPTION_TIME ).total_seconds()
                minutes_since_eruption.append( seconds_since_eruption / 60 )

                seconds_since_previous = ( current_timestamp - previous_timestamp).total_seconds()
                pd_hPa_minute = ( current_pressure_hPa - previous_pressure_hPa ) * 60.0 / seconds_since_previous
                _pressure_diff_hPa_minute.append( pd_hPa_minute )

            previous_timestamp    = current_timestamp
            previous_pressure_hPa = current_pressure_hPa

        return _pressure_diff_hPa_minute, minutes_since_eruption

    def pressure_diff_hPa_minute_with_timestamp( self, interpolate = None ):
        _pressure_diff_hPa_minute, _minutes_since_eruption = \
            self.pressure_diff_hPa_minute( interpolate )
        
        return \
            _pressure_diff_hPa_minute, \
            [ ERUPTION_TIME + timedelta( minutes = d_minute )
              for d_minute in _minutes_since_eruption ]
        

class BarometricPressureRecord:
    '''One record of barometric pressure data from one monitoring station.'''

    TABLE_NAME = 'barometric_pressure_data'
    
    def __init__(
        self, 
        station,
        timestamp,
        pressure_hPa
    ):
        self.station      = station
        self.timestamp    = timestamp
        self.pressure_hPa = pressure_hPa

    @classmethod
    def create_table( cls, database ):
        '''
        Create the table of barometric pressure data in the given SQLite3 database
        if not exists.
        '''

        if database.table_exists( cls.TABLE_NAME ):
            return
        
        tqdm.write( f'Creating table "{cls.TABLE_NAME}" ...' )
        
        sql_statement = f'''
CREATE TABLE IF NOT EXISTS {cls.TABLE_NAME} (
    data_source   TEXT NOT NULL,
    station_id    TEXT NOT NULL,
    timestamp     TIMESTAMP NOT NULL,
    pressure_hPa  REAL NOT NULL,

    FOREIGN KEY (data_source) REFERENCES {BarometricPressureMonitoringStation.TABLE_NAME}(data_source),
    FOREIGN KEY (station_id)  REFERENCES {BarometricPressureMonitoringStation.TABLE_NAME}(identifier),
    PRIMARY KEY (data_source, station_id, timestamp)
);
'''
        database.execute( sql_statement )

    @classmethod
    def insert_records( cls, database, records ):
        database.executemany( f'''
INSERT OR IGNORE INTO {cls.TABLE_NAME}
(data_source, station_id, timestamp, pressure_hPa)
VALUES
(?, ?, ?, ?)
    ''', [ ( record.station.data_source,
             record.station.identifier,
             record.timestamp.astimezone( timezone.utc ),
             record.pressure_hPa) for record in records ] )

