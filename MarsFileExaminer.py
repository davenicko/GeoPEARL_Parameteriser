# -*- coding: utf-8 -*-
import pandas as pd

class climate_mars:
    '''
    Class to handle the MARS database output in order to extract relevant
    information (lat, long of stations etc.)
    '''
    
    def __init__(self, file_path):
        '''
        Initialise the class. Read the met file and change the day to a
        datetime format.

        file_path:  path to the MARS met file to import
        '''
        self.clim_data = pd.read_csv(file_path, na_values=0, delimiter=';')
        self.clim_data['DAY'] = pd.to_datetime(self.clim_data['DAY'],
                                               format='%Y%m%d')
        self.coordinate_list = []
        print('Initialisaton complete...')

    def __repr__(self):
        '''
        return the repr of the dataframe of the imported data
        '''
        return repr(self.clim_data)

    def get_lat_long(self):
        '''
        Retrive a list of latitude and longitude coordinates and return a
        dictionary containing the results
        '''
        return self.clim_data.groupby(['GRID_NO', 'LATITUDE', 'LONGITUDE'])\
               .size().reset_index().rename(columns={0:'count'})

    def output_stations(self, output_file):
        '''
        Output the list of stations to a csv file
        '''
        self.coords = self.get_lat_long()
        self.coords.to_csv(output_file, encoding='utf-8')

    def clip_data_to_region(self, clipped_locations):
        '''
        Clip the data to include only the locations in clipped_locations
        '''
        stations = pd.read_csv(clipped_locations, na_values=0)
        mask = self.clim_data['GRID_NO'].isin(stations['GRID_NO'].values)
        self.clipped_clim_data = self.clim_data.loc[mask]

    def clip_data_to_date(self, start_date, end_date):
        '''
        Return a dataframe containing the data between only two dates
        '''
        mask = (self.clim_data['DAY']>=start_date) &\
               (self.clim_data['DAY']<=end_date)
        self.clipped_clim_data_date = self.clipped_clim_data.loc[mask]

    def export_clipped_data(self, output_file):
        '''
        Export the clipped data to csv
        '''
        self.clipped_clim_data_date.to_csv(output_file, encoding='utf8')

    def rearrange_dataframe(self):
        '''
        Rearrange the data frame ready for export
        '''
        self.clipped_clim_data_date['day'] = self.clipped_clim_data_date['DAY'].dt.day
        self.clipped_clim_data_date['month'] = self.clipped_clim_data_date['DAY'].dt.month
        self.clipped_clim_data_date['year'] = self.clipped_clim_data_date['DAY'].dt.year
        self.clipped_clim_data_date['rad'] = -99
        self.clipped_clim_data_date['hum'] = -99
        self.clipped_clim_data_date['wind'] = -99
        self.clipped_clim_data_date['tmin'] = self.clipped_clim_data_date['TEMPERATURE_AVG']
        self.clipped_clim_data_date['tmax'] = self.clipped_clim_data_date['TEMPERATURE_AVG']
        columns = ['GRID_NO', 'day', 'month', 'year', 'rad', 'tmin', 'tmax',
                   'hum', 'wind', 'PRECIPITATION', 'ET0']
        return self.clipped_clim_data_date.loc[:,columns]

    def metfile_creator(self, output_folder):
        '''
        Function to create a metfile for each station in the input file and save tothe given
        output file.
        '''

        Header = ('*--------------------------------------------------------------\n'
                  '* Meteofiles for GeoPEARL\n'
                  '* Version 2 - ETREF multiplied by 0.8 to obtain Makkink\n'
                  '*\n'
                  '* Station DD MM YYYY RAD   Tmin Tmax HUM WIND RAIN ETref\n'
                  '*         nr nr nr   kJ/m2 C    C    kPa m/s  mm   mm\n'
                  '*--------------------------------------------------------------\n')
        
        for station in self.clipped_clim_data_date['GRID_NO'].unique():
            mask = self.clipped_clim_data_date['GRID_NO']==station
            climate_output = self.rearrange_dataframe().loc[mask]
            try:
                climate_output.iloc[0:3, 0:3]
                climate_output['PRECIPITATION'] = climate_output['PRECIPITATION'].fillna(0)
                if output_folder[-1] == '\\':
                    output_file = output_folder + str(station) + ".met"
                else:
                    output_file = output_folder + '\\' + str(station) + ".met"
                
                with open(output_file, 'w') as file_header:
                    file_header.write(Header)
    
                climate_output.to_csv(output_file, mode= 'a', header=None,
                                      sep=' ', encoding='utf8', index=False)
            except:
                print(station + " Failed")
            
