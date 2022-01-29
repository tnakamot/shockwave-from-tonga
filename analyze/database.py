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
             record.timestamp,
             record.pressure_hPa) for record in records ] )

