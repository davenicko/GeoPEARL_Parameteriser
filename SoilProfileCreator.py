# -*- coding: utf-8 -*-
import os
import math
import csv
from collections import namedtuple
import logging
import datetime
import pandas as pd
import numpy as np
from timeit import default_timer as timer


class Soil_profile_creator:
    '''
    Class to take in data from SPADE to link to soil typologiacal units (STUs)
    and generate the soil profiles for use in GeoPEARL using pedo transfer
    functions (PTFs) where necessary

    Should hold a dataframe with the columns in the same order as the GeoPEARL
    .sol file

    This can then be exported with the function sol_export to generate the .sol
    file
    '''
    def __init__(self,
                 stu_to_prof_name_file,
                 prof_name_to_int_file,
                 spade_hor_file,
                 toth_asc_dir,
                 hypres_asc_dir,
                 output_folder,
                 mode,
                 ):
        '''
        Initialise the data

        stu_to_prof_name_file: csv file linking the regions STUs to profile names
        prof_name_to_int_file: csv file linking profile names to integer repr
                               N.B. some prof names may be NULL
        spade_hor_file     : File containing the soil database in excel format
        toth_asc_dir         : directory containing MVG paramter rasters
        mode                 : Which mode to use, determines the pedo-
                               transfer functions used
                               Options: HYPRES, or Toth

        returns: nothing, but sets up the object to output a .sol file
        '''
        logging.basicConfig(filename='GeoPEARLparameteriser.log', level=logging.WARNING)
        self.border = ('\n========================================'
                       '========================================\n')
        self.stu_prof = pd.read_csv(stu_to_prof_name_file)
        prof_int_df = pd.read_csv(prof_name_to_int_file)
        self.prof_int = dict(zip(prof_int_df.fid, prof_int_df.Links_prof_stu_PROF_NUM))
        self.spade_hor = pd.read_excel(spade_hor_file)
        self.toth_asc_dir = toth_asc_dir
        self.hypres_asc_dir = hypres_asc_dir
        self.output_folder = output_folder
        self.mode = mode
        self.func_time = 0
        self.scenario_hashes = []
        self.scaling_timer = {}
        self.plo_dict_to_export = {}
        self.scaling_timer_start = timer()
        self.scen_met_list = self.rast_extract(hypres_asc_dir,
                                               'NearestMetStation',
                                               skip=6)

        '''
        =======================================================================
        Checking variables to delete
        =======================================================================
        '''
        self.param_n_check = []

        #Default values
        self.res_wat_cont = 0
        #If needed for soil risk assessment thie dispersion length can be set
        #to 2.5cm
        self.Disp_length = 0.05 #See Tiktak 2012
        self.Sesqui_ox_cont = 200

        #list of data in order of columns
        self.data_list = [
            'PROF_NUM',
            'HOR_NUM',
            'HOR_THICK',
            'NUM_NUMER_COMP',#See Tiktak 2012
            'SAND_FRAC',
            'SILT_FRAC',
            'CLAY_FRAC',
            'OM',
            'PH_KCL',
            'SAT_SOL_WAT_CONT', #Theta s
            'RES_WAT_CONT', #Theta r
            'PARAM_ALPHA_DRY',
            'PARAM_ALPHA_WET',#Same as above
            'PARAM_N',
            'SAT_HYD_CON', #Ks
            'PHY_SAT_HYD_CON',#Same as above
            'PARAM_L',
            'DISP_LEN',
            'SESQ_OX_CONT',
            ]

    def __repr__(self):
        '''
        repr for the data. Return the repr for the spade file.
        '''
        return repr(self.spade_hor)

    def sol_creator(self):
        '''
        Uses the spade file to output a file with the relevant soil
        properties including those transformed using the relevan PTFs

        This will output for an entire country in Europe which can then be
        clipped to the region of interest
        '''
        #Generate the ptfs according to the wanted method
        if self.mode == 'HYPRES':
            start = timer()
            self.HYPRES_PTF()
            end = timer()
            start_plo = timer()
            self.hypres_plot_creator()
            end_plo = timer()

            logging.warning('self.HYPRES_PTF() took %3.5f seconds to run' % 
                            (end-start))
            logging.warning('self.hypres_plot_creator() took %3.5f seconds to '
                            'run' % (end_plo-start_plo))
            logging.warning(self.border + 'Run ends' + self.border)
        else:
            start = timer()
            self.Toth_PTF()
            end = timer()

            logging.warning('self.Toth_PTF() took %3.5f seconds to run' % (end-start))
            logging.warning('toth_spade_linker part took %3.5f seconds to run' % self.func_time)
            logging.warning(self.border + 'Run ends' + self.border)
        
    def HYPRES_PTF(self):
        '''
        Uses the spade file to output a file with the relevant soil
        properties including those transformed using the HYPRES PTFs

        '''

        def ptf_func(row, func):
            '''
            PTF to determine saturated soil water content
            '''
            C = row['TEXT_2']
            S = row['TEXT_20'] + row['TEXT_50']
            OM = row['OM']
            D = row['BD']
            if row['DEPTHSTART'] < 30: #assuming topsoil is upper 30cm
                topsoil = 1
            else:
                topsoil = 0
            return func(C, S, OM, D, topsoil)

        def theta_s(C, S, OM, D, topsoil):
            '''
            Determine the theta_s parameter using the HYPRES PTFs
            '''
            if C > 0 and S > 0 and OM > 0 and D > 0:
                ptf_res = 0.7919 + 0.001691*C - 0.29619*D - 0.000001491*S*S\
                        + 0.0000821*OM*OM + 0.02427*(1/C) + 0.01113*(1/S)\
                        + 0.01472*math.log(S) - 0.0000733*OM*C - 0.000619*D*C\
                        - 0.001183*D*OM - 0.0001664*topsoil*S
                return ptf_res
            else:
                return np.NAN

        def param_alpha(C, S, OM, D, topsoil):
            '''
            Determine the alpha parameter using the HYPRES PTFs
            '''
            if C > 0 and S > 0 and OM > 0 and D > 0:
                ptf_res = -14.96 + 0.03135*C + 0.0351*S + 0.646*OM +15.29*D\
                        - 0.192*topsoil - 4.671*D*D - 0.000781*C*C - 0.00687*OM*OM\
                        + 0.0449*(1/OM) + 0.0663*math.log(S) + 0.1482*math.log(OM)\
                        - 0.04546*D*S - 0.4852*D*OM +0.00673*topsoil*C
                #return the parameter with the transformation reversed
                return math.exp(ptf_res)
            else:
                return np.NAN

        def param_n(C, S, OM, D, topsoil):
            '''
            Determine the n parameter using the HYPRES PTFs
            '''
            if C > 0 and S > 0 and OM > 0 and D > 0:
                ptf_res =-25.23 - 0.02195*C + 0.0074*S - 0.1940*OM + 45.5*D\
                        - 7.24*D*D + 0.0003658*C*C + 0.002885*OM*OM -12.81*(1/D)\
                        - 0.1524*(1/S) - 0.01958*(1/OM) - 0.2876*math.log(S)\
                        - 0.0709*math.log(OM) - 44.6*math.log(D) - 0.02264*D*C\
                        + 0.0896*D*OM + 0.00718*topsoil*C
                #return the parameter with the transformation reversed
                return math.exp(ptf_res)+1
            else:
                return np.NAN

        def param_l(C, S, OM, D, topsoil):
            '''
            Determine the L parameter using the HYPRES PTFs
            '''
            if C > 0 and S > 0 and OM > 0 and D > 0:
                ptf_res = 0.0202 + 0.0006193*C*C - 0.001136*OM*OM\
                        - 0.2316*math.log(OM) - 0.03544*D*C + 0.00283*D*S\
                        + 0.0488*D*OM
                #return the parameter with the transformation reversed
                return 10*((math.exp(ptf_res)-1)/(math.exp(ptf_res)+1))
            else:
                return np.NAN

        def Ks(C, S, OM, D, topsoil):
            '''
            Determine the Ks parameter using the HYPRES PTFs
            '''
            if C > 0 and S > 0 and OM > 0 and D > 0:
                ptf_res = 7.755 + 0.0352*S + 0.93*topsoil - 0.967*D*D\
                        - 0.000484*C*C - 0.000322*S*S + 0.001*(1/S) - 0.0748*(1/OM)\
                        - 0.643*math.log(S) - 0.01398*D*C - 0.1673*D*OM\
                        + 0.02986*topsoil*C - 0.03305*topsoil*S
                #return the parameter with the transformation reversed
                return math.exp(ptf_res)/100
            else:
                return np.NAN

        self.PTF_df = pd.DataFrame()
        self.PTF_df['SAT_SOL_WAT_CONT'] =\
                self.spade_hor.apply(lambda row: ptf_func(row, theta_s), axis = 1)
        self.PTF_df['PARAM_ALPHA_DRY'] =\
                self.spade_hor.apply(lambda row: ptf_func(row, param_alpha), axis = 1)
        self.PTF_df['PARAM_N'] =\
                self.spade_hor.apply(lambda row: ptf_func(row, param_n), axis = 1)
        self.PTF_df['PARAM_L'] =\
                self.spade_hor.apply(lambda row: ptf_func(row, param_l), axis = 1)
        self.PTF_df['SAT_HYD_CON'] =\
                self.spade_hor.apply(lambda row: ptf_func(row, Ks), axis = 1)

        temp_sol_df = pd.DataFrame()
        temp_sol_df['PROF_NUM'] = self.spade_hor['PROF_NUM'].str[2:]
        temp_sol_df['HOR_NUM'] = self.spade_hor['HOR_NUM']
        HOR_THICK = self.spade_hor['DEPTHEND'] - self.spade_hor['DEPTHSTART']
        temp_sol_df['HOR_THICK'] = HOR_THICK / 100.0
        temp_sol_df['NUM_NUMER_COMP'] = \
                self.spade_hor.apply(lambda row: self.num_numeric(row), axis=1)
        temp_sol_df['SAND_FRAC'] = self.spade_hor['TEXT_2000'] +\
                                   self.spade_hor['TEXT_200']
        temp_sol_df['SILT_FRAC'] = self.spade_hor['TEXT_20'] +\
                                   self.spade_hor['TEXT_50']
        temp_sol_df['CLAY_FRAC'] = self.spade_hor['TEXT_2']
        temp_sol_df['OM'] = self.spade_hor['OM']
        temp_sol_df['PH_KCL'] = self.spade_hor['PH']
        temp_sol_df['SAT_SOL_WAT_CONT'] = self.PTF_df['SAT_SOL_WAT_CONT']
        temp_sol_df['RES_WAT_CONT'] = self.res_wat_cont
        temp_sol_df['PARAM_ALPHA_DRY'] = self.PTF_df['PARAM_ALPHA_DRY']
        temp_sol_df['PARAM_ALPHA_WET'] = self.PTF_df['PARAM_ALPHA_DRY']#same as above
        temp_sol_df['PARAM_N'] = self.PTF_df['PARAM_N']
        temp_sol_df['SAT_HYD_CON'] = self.PTF_df['SAT_HYD_CON']
        temp_sol_df['PHY_SAT_HYD_CON'] = self.PTF_df['SAT_HYD_CON']#same as above
        temp_sol_df['PARAM_L'] = self.PTF_df['PARAM_L']
        temp_sol_df['DISP_LEN'] = self.Disp_length
        temp_sol_df['SESQ_OX_CONT'] = self.Sesqui_ox_cont

        #change empty parameters to NaN:
        temp_sol_df = temp_sol_df.replace(np.nan, -999)
        temp_sol_df = temp_sol_df[temp_sol_df.astype(int) > -900]
        
        #select only the soil profiles that are in the study area
        #First get a list of soil profiles (as integers)
        sol_prof_list = []
        for soil in temp_sol_df.PROF_NUM:
            prof_num = self.get_key(self.prof_int, 'FR' + soil)
            if prof_num and prof_num not in sol_prof_list:
                sol_prof_list.append(prof_num)

        #Now use that list to filter out the profiles that are in the study area
        self.sol_df = pd.DataFrame()
        for prof_num in sol_prof_list:
            temp_line = temp_sol_df.loc[temp_sol_df.PROF_NUM == self.prof_int[prof_num][2:]]
            temp_line.PROF_NUM = prof_num
            self.sol_df = self.sol_df.append(temp_line)


    def hypres_plot_creator(self):
        '''
        Method to create plots from the scenario and met rasters
        '''
        hypres_files = ['NearestMetStation', 'ProfileNameRaster']
        self.hypres_plot_raster = []
        asc_dict = {}
        self.plot_dict = {}
        self.plot_area_dict = {}
        for file in hypres_files:
            asc_dict[file] = self.rast_extract(self.hypres_asc_dir, file, skip=6)
        met_file_gen = self.param_gen(asc_dict['NearestMetStation'])
        #For each line in the rasters:
        #For each cell in the line, if the cell value == 0 there is no soil
        #profile so move on. If there is a soil profile check the associated
        #met file and if this combination is unique, add to a list of plots.
        #if not don't add to the list, but either way add the plot number to
        #a list
        #at the end of each line add the list of plots to a list which will
        #be exported as an asc
        plot_num = 1
        for line in asc_dict['ProfileNameRaster']:
            plot_raster_line = []
            for cell in line:
                if self.mode =='HYPRES':
                    met = next(met_file_gen)
                    if cell:
                        if (cell, met) not in self.plot_dict.values():
                            self.plot_dict[plot_num] = (cell, met)
                            self.hypres_plot_raster.append(plot_num)
                            self.plot_area_dict[plot_num] = 1
                            plot_num += 1
                        else:
                            self.plot_area_dict[self.get_key(self.plot_dict, (cell, met))] += 1
                            self.hypres_plot_raster.append(self.get_key(self.plot_dict, (cell, met)))
                    else:
                        self.hypres_plot_raster.append('-')

                else:
                    #create a temp subset of the sol data for the soil profile:
                    met = next(met_file_gen)
                    if cell:
                        try:
                            prof = int(self.prof_int[int(cell)][2:])
                        except KeyError:
                            prof = -99
                        temp_sol_df = self.sol_df[self.sol_df.PROF_NUM == prof]

                        if len(temp_sol_df) and temp_sol_df.SAND_FRAC.iloc[0]:
                            if (cell, met) not in self.plot_dict.values():
                                self.plot_dict[plot_num] = (cell, met)
                                plot_raster_line.append(plot_num)
                                self.plot_area_dict[plot_num] = 1
                                plot_num += 1
                            else:
                                self.plot_area_dict[self.get_key(self.plot_dict, (cell, met))] += 1
                                plot_raster_line.append(self.get_key(self.plot_dict, (cell, met)))
                        else:
                            plot_raster_line.append('-')                
                    else:
                        plot_raster_line.append('-')                



    def plo_creator_hypres(self):
        '''
        Method to create the plo dataframe ready for export

        self.plot_dict

        is the list needed to cycle through
        '''
        crop_type = 1
        cor_precip = 1
        cor_temp = 0
        cor_et = 1
        irr_sw = 0
        max_pond_dep = 0.01
        air_boun_th = 0.01
        vuln_rank = 1 #For now give all scenarios the same vulnerability rank
        plo_toexp = {}
        plot = 1
        columns = ['plot', 'area', 'meteo', 'soil', 'crop', 'cor_precip',
                   'cor_temp', 'cor_et', 'irr_sw', 'max_pond_dep',
                   'air_boun_th', 'vuln_rank', 'text_class']
        self.toth_plot_raster = [] 
        for col in columns:
            plo_toexp[col] = []

        try:
            for prof in self.plot_dict.keys():
                #calculate the soil profile / met combinations (plot == scenario)
                om = self.sol_df.loc[self.sol_df.PROF_NUM == prof].OM.mean()*100
                sand = self.sol_df.loc[self.sol_df.PROF_NUM == prof].SAND_FRAC.mean()*100
                clay = self.sol_df.loc[self.sol_df.PROF_NUM == prof].CLAY_FRAC.mean()*100
                plo_toexp['plot'].append(plot)
                plo_toexp['area'].append(10)
                plo_toexp['meteo'].append(self.plot_dict[prof][1])
                plo_toexp['soil'].append(self.plot_dict[prof][0])
                plo_toexp['crop'].append(crop_type)
                plo_toexp['cor_precip'].append(cor_precip)
                plo_toexp['cor_temp'].append(cor_temp)
                plo_toexp['cor_et'].append(cor_et)
                plo_toexp['irr_sw'].append(irr_sw)
                plo_toexp['max_pond_dep'].append(max_pond_dep)
                plo_toexp['air_boun_th'].append(air_boun_th)
                plo_toexp['vuln_rank'].append(vuln_rank)
                plo_toexp['text_class'].append(self.text_calc(sand, clay, om))
                self.toth_plot_raster.append(str(plot))
                plot += 1

        except Exception:
            print('Error - was sol_creator() executed?')

        self.plo_dict_to_export = pd.DataFrame(plo_toexp, columns=columns)

    def plo_creator_toth(self):
        '''
        Method to create the plo dataframe ready for export
        
        self.scen_raster_list
        
        and
        
        self.scen_met_list
        
        are the lists needed to cycle through
        '''
        crop_type = 1
        cor_precip = 1
        cor_temp = 0
        cor_et = 1
        irr_sw = 0
        max_pond_dep = 0.01
        air_boun_th = 0.01
        vuln_rank = 1 #For now give all scenarios the same vulnerability rank
        plot_dict = {}
        plot = 1
        soil_prof_gen = self.param_gen(self.scen_raster_list)
        columns = ['plot', 'area', 'meteo', 'soil', 'crop', 'cor_precip',
                   'cor_temp', 'cor_et', 'irr_sw', 'max_pond_dep',
                   'air_boun_th', 'vuln_rank', 'text_class']
        plot_dict_check = {}
        plot_areas = {}
        
        #a list to export as a raster containing the plot numbers
        self.toth_plot_raster = [] 
        for col in columns:
            plot_dict[col] = []
        met_gen = self.param_gen(self.scen_met_list)

        try:
            for prof in soil_prof_gen:
                met = next(met_gen)
                if prof:
                    #calculate the soil profile / met combinations (plot == scenario)
                    if (prof, met) not in list(plot_dict_check.keys()):
                        om = self.sol_df.loc[self.sol_df.PROF_NUM ==1].OM.mean()*100
                        sand = self.sol_df.loc[self.sol_df.PROF_NUM ==1].SAND_FRAC.mean()*100
                        clay = self.sol_df.loc[self.sol_df.PROF_NUM ==1].CLAY_FRAC.mean()*100
                        plot_dict['plot'].append(plot)
                        plot_dict['area'].append(10)
                        plot_dict['meteo'].append(met)
                        plot_dict['soil'].append(prof)
                        plot_dict['crop'].append(crop_type)
                        plot_dict['cor_precip'].append(cor_precip)
                        plot_dict['cor_temp'].append(cor_temp)
                        plot_dict['cor_et'].append(cor_et)
                        plot_dict['irr_sw'].append(irr_sw)
                        plot_dict['max_pond_dep'].append(max_pond_dep)
                        plot_dict['air_boun_th'].append(air_boun_th)
                        plot_dict['vuln_rank'].append(vuln_rank)
                        plot_dict['text_class'].append(self.text_calc(sand, clay, om))
                        plot_dict_check[(prof, met)] = plot
                        plot_areas[plot] = 1
                        print('Plot %s stored' % plot)
                        self.toth_plot_raster.append(str(plot))
                        plot += 1

                    else:
                        self.toth_plot_raster.append(plot_dict_check[(prof, met)])
                        plot_areas[plot_dict_check[(prof, met)]] += 1
                else:
                    self.toth_plot_raster.append('-')


        except Exception:
            print('Error - was sol_creator() executed?')

        plot_dict['area'] = pd.Series(plot_areas.values())

        #Make sure that the plots are in the study area
        self.plo_dict_to_export = pd.DataFrame(plot_dict, columns=columns)
  
    def text_calc(self, sand, clay, om):
        '''
        Test for soil type (sand, clay or loam) according to the UK texture
        triangle
        1 = sand, 2 = clay, 3 = loam, 4 = peat
        '''
        #Test if the soil is Peat. From Natural England Technical Information
        #Note TIN037
        if om > 25:
            return 4
        if (om - 0.1 * clay) > 20:
            return 4
        
        #If not Peat, test for soil texture
        if sand - clay > 70:
            return 1
        if clay >= 35:
            return 2
        if 50 > sand > 45 and clay + sand >= 80:
            return 2
        if sand >= 50 and clay >= 30:
            return 2
        return 3

    def num_numeric(self, row):
        '''
        According to Tiktak 2012 the following numerical soil components
        are suitable for the given depths

        If used for soil risk assessment the thickness of the numerical
        layers should be set to smaller values for the top 1cm (see EFSA,
        2010)
        '''
        HOR_THICK = row['DEPTHEND'] - row['DEPTHSTART']
        if row['DEPTHEND'] <= 30:
            return HOR_THICK / 1.0
        elif row['DEPTHEND'] <= 100:
            return HOR_THICK / 2.5
        else:
            return HOR_THICK / 5.0

    def value_check(self, value, limit):
        '''
        Check a value is above a prescribed limit
        '''
        if value > limit:
            return value
        return np.nan

    def Toth_PTF(self):
        '''
        Uses the toth ptf files to output a the relevant soil
        properties including those transformed using the Toth PTFs

        The toth files are named with a letter and a number
        Letter key:
            a = alpha               \n
            k = K0                  \n
            L = L/Lambda            \n
            n = n                   \n
            t = ThetaS              \n
        The numbers denote depth:
           1 = sl1 = 0 cm           \n
           2 = sl2 = 5 cm           \n
           3 = sl3 = 15 cm          \n
           4 = sl4 = 30 cm          \n
           5 = sl5 = 60 cm          \n
           6 = sl6 = 100 cm         \n
           7 = sl7 = 200 cm         \n
        Each file contains the study area in raster format projected to the
        Lambert_Azimuthal_Equal_Area projection. The area is in square format
        and therefore there are values to represent areas outside the clipped
        raster (-2147483647 for most, extreme edges are 0). Values are x10,000
        so need to be divided to a float
        '''
        msg = self.border + 'Toth run began ' + str(datetime.datetime.now()) + self.border
        logging.warning(msg)
        toth_files = [
                'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7',
                'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k7',
                'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7',
                'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7',
                't1', 't2', 't3', 't4', 't5', 't6', 't7',
                ]

        self.raster_dict = {}
        for file_name in toth_files:
            self.raster_dict[file_name] = self.rast_extract(self.toth_asc_dir,
                                                            file_name, 
                                                            dtype=float)/10000.0
        #Add in the ProfRast separately
        self.raster_dict['ProfRast'] = self.rast_extract(self.toth_asc_dir,
                                                         'ProfRast')
        self.raster_dict['NearestMetStation'] = self.rast_extract(self.toth_asc_dir, 
                                                                  'NearestMetStation',
                                                                  skip=6)

        #For printing the progress of the run we need to know the size of the
        #study area in number of square kilometers:
        self.area_size = self.raster_dict['ProfRast'].size

        #Append these to the list of files
        toth_files.append('ProfRast')
        toth_files.append('NearestMetStation')

        #Create a dict with an entry for each unique profile. This dict
        #contains all the parameters to enable running of the groundwater model
        #and is numbered sequentially from 0 onwards
        profile_dict = {}


        #for each pixel, take all the parameters needed and store them in a
        #named tuple. Store these named tuples in a dict. If the named tupple
        #is not in the list then add the tuple to the profile_dict.
        #Add the number of the pixel to a dict with the hash as a key so we
        #know which cell contains which profile
        toth_params = namedtuple('toth_params',
                                 (
                                  'HOR_NUM '
                                  'HOR_THICK '
                                  'NUM_NUMER_COMP '#See Tiktak 2012
                                  'SAND_FRAC '
                                  'SILT_FRAC '
                                  'CLAY_FRAC '
                                  'OM '
                                  'PH_KCL '
                                  'SAT_SOL_WAT_CONT ' # = t1-7
                                  'RES_WAT_CONT ' # = 0
                                  'PARAM_ALPHA_DRY ' # = a1-7
                                  'PARAM_ALPHA_WET '# = a1-7
                                  'PARAM_N ' # = n1-7
                                  'SAT_HYD_CON ' # = k1-7
                                  'PHY_SAT_HYD_CON '# = k1-7
                                  'PARAM_L ' # = L1-7
                                  'DISP_LEN ' # = 0.05
                                  'SESQ_OX_CONT ' # = 200
                                  #'NearestMetStation '
                                  )
                                 )

        #Set up a tuple containing all zeros for the out of scope areas
        self.zero_params = toth_params(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0,)
        #store the generators for each parameter in a dictionary so they can be
        #easily accessed
        gen_store = {}
        for file in toth_files:
            gen_store[file] = self.param_gen(self.raster_dict[file])

        depth_list = ['1', '2', '3', '4', '5', '6', '7']
        parameter_list = ['t', 'a', 'n', 'k', 'L']
        thicknesses = [5, 10, 15, 30, 40, 100, 100]
        num_num = [5, 10, 15, 12, 16, 20, 20]
        thick_dict = dict(zip(depth_list, thicknesses))
        num_num_dict = dict(zip(depth_list, num_num))
        cell_no = 0
        scenario = 0
        scenario_list = []
        
        #Lists containing a raster with the scenarios or met station for each
        #raster square
        self.scen_raster_list = []
        self.scen_met_list = []
        #List for each line of these rasters
        scen_raster_line =[]
        scen_met_line = []

        self.test_counter = 0

        for row in self.raster_dict['ProfRast']:
            for cell in row:
                profiles = []
                nearest_met_station = next(gen_store['NearestMetStation'])
                scen_met_line.append(nearest_met_station)

                for depth in depth_list:
                    hydr_props = {}
                    for param in parameter_list:
                        hydr_props[param] = (next(gen_store[param+depth]))
                    if cell:
                        #Sometimes the cell may inlcude a Null value in which case
                        #this will fail. The try block catches this but also may
                        #mask other errors.
                        try:
                            spade_props = self.toth_spade_linker(cell, int(depth))

                            soil_prof = toth_params(
                                    int(depth),
                                    thick_dict[depth],
                                    num_num_dict[depth],
                                    spade_props.SAND_FRAC,
                                    spade_props.SILT_FRAC,
                                    spade_props.CLAY_FRAC,
                                    spade_props.OM,
                                    spade_props.PH_KCL,
                                    hydr_props['t'],
                                    0,
                                    hydr_props['a'],
                                    hydr_props['a'],
                                    hydr_props['n'],
                                    hydr_props['k'],
                                    hydr_props['k'],
                                    hydr_props['L'],
                                    0.05,
                                    200,
                                    )

                        except Exception:
                            logging.warning('Run with soil ID ' + str(cell)\
                                            + ' Failed')
                            soil_prof = toth_params(int(depth), 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0,)
                    else:
                        soil_prof = toth_params(int(depth), 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0,)

                    profiles.append(soil_prof)

                if soil_prof.PARAM_N not in self.param_n_check:
                    self.param_n_check.append(soil_prof.PARAM_N)
                soil_prof_hash = hash(tuple(profiles))

                #Check to see if the first depth has a valid value for N. It is
                #certain that as N has a value > 1, then valid values for N in
                #this raster are < 10,000. Therefore if values are < 1, i.e.
                #0 or negative, then ignore them
                if profiles[0].PARAM_N > 1:
                    if soil_prof_hash not in self.scenario_hashes:
                        scen_raster_line.append(scenario)
                        self.scenario_hashes.append(soil_prof_hash)
                        profile_dict[scenario] = tuple(profiles)
                        for i in range(7):
                            scenario_list.append(str(scenario))
                        scenario += 1
                    else:
                       scen_raster_line.append(self.get_key(profile_dict,
                                                            tuple(profiles)))
                else:
                    #If the value for N is not valid, pad with an empty tuple
                    #for each depth
                    soil_prof = toth_params(int(depth), 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0,)
                    profiles = []
                    for i in range(7):
                        profiles.append(soil_prof)
                    scen_raster_line.append(0)
                cell_no += 1

#==============================================================================
                #Progress code - so we know how far through we are (this takes
                #some time)
                if cell_no % 100 == 0:
                    percent = cell_no/self.area_size*100
                    print('Progress = %3.2f percent' % percent)

#==============================================================================

            self.scen_raster_list.append(scen_raster_line)
            scen_raster_line = []
            self.scen_met_list.append(scen_met_line)
            scen_met_line = []

        #Add the data to a df for export. We have all paramaters in a tuple
        #containing one tuple for each cell that contains all 7 depths as
        #prescribed in the Toth PTFs

        self.toth_df_creator(profile_dict, scenario_list)

    def toth_df_creator(self, profile_dict, scenario_list):
        columns = ['HOR_NUM',
                   'HOR_THICK',
                   'NUM_NUMER_COMP',
                   'SAND_FRAC',
                   'SILT_FRAC',
                   'CLAY_FRAC',
                   'OM',
                   'PH_KCL',
                   'SAT_SOL_WAT_CONT',
                   'RES_WAT_CONT',
                   'PARAM_ALPHA_DRY',
                   'PARAM_ALPHA_WET',
                   'PARAM_N',
                   'SAT_HYD_CON',
                   'PHY_SAT_HYD_CON',
                   'PARAM_L',
                   'DISP_LEN',
                   'SESQ_OX_CONT',
                   ]
        self.sol_df = pd.DataFrame(columns=columns)
        for scenario in profile_dict.keys():
            for depth in profile_dict[scenario]:
                temp_dict = dict(zip(columns, depth))
                self.sol_df = self.sol_df.append(temp_dict, ignore_index=True)
            prof_num_series = pd.Series(scenario_list)
            self.sol_df['PROF_NUM'] = prof_num_series


    def get_key(self, sample_dict, val):
        '''
        return a key from a dict given a value
        '''
        for key, value in sample_dict.items():
            if val == value:
                return key

    def toth_spade_linker(self, prof_num, hor_num):
        '''
        =======================================================================
        This is the bottleneck. About 75% of the running time is spent in this
        function when using a small sample of data. If a larger sample is used
        then this takes exponentiallly longer as it searched for clashes in the
        scenario data structure. This takes about 35 minutes for the Pay de la
        Loire region.

        Ideas for speeding up in the future:

            https://towardsdatascience.com/how-to-make-your-pandas-loop-71-803-times-faster-805030df4f06
        =======================================================================
        Link the soil profile number of the raster with the SPADE row that
        coincides with it. The data is then exracted and returned from this row.
        Note that the depths in the 3D hydraulic database of Europe are
        consistent, but the depths of the SPADE are not therefore this function
        finds the horizon that contains most of the Toth soil column and uses
        that data. If there is a tie it defaults to the lower depth.

        prof_num: the number of the profile to link with spade
        hor_num: the horizon number (1-7)

        Return a tuple containing: profile number,
                                   horizon number,
                                   horizon thickness,
                                   number of numerical compartments,
                                   sand fraction,
                                   silt fraction,
                                   clay fraction,
                                   organic matter content,
                                   pH in KCl
        '''
        horizon_depths = {1: 0, 2: 5, 3: 15, 4: 30,
                          5: 60, 6: 100, 7: 200, 8: 1000}
        profile_name = self.prof_int[prof_num]

        #Timing for future optimisation based on a ~4 second run time
        ######## 0.21 - 0.28s ########
        profile_df = pd.DataFrame(
                self.spade_hor.loc[self.spade_hor.PROF_NUM == profile_name])

        depth_start = horizon_depths[hor_num]
        depth_end = horizon_depths[hor_num + 1]
        ######## end ########

        func_start = timer()
        #calculate the percentage the toth depth is in the spade depth
        ######## 0.44 - 0.45s ######## - could use pandas vectorisation?
        profile_df['perc_col_in_hor'] = profile_df.apply(
                lambda row: self.hor_perc_calc(depth_start, depth_end, row),
                axis=1)
        func_end = timer()
        self.func_time += func_end - func_start
        ######## end ########

        ######## 0.4s ########
        #3.3m or 12% of a 28 min run
        if profile_df.perc_col_in_hor.any():
            #Minimize the dataframe to only those rows with the highest percent
            #and if there is a tie we will take the max depth
            profile_df = pd.DataFrame(
                    profile_df.loc[profile_df.perc_col_in_hor ==
                                   profile_df.perc_col_in_hor.max()])
        if len(profile_df.index) > 1:
            profile_df = pd.DataFrame(profile_df.loc[profile_df.HOR_NUM ==
                                                     profile_df.HOR_NUM.max()])
        ######## end ########


        spade_data = namedtuple('spade_data',(
                                'PROF_NUM '
                                'HOR_NUM '
                                'HOR_THICK '
                                'NUM_NUMER_COMP '
                                'SAND_FRAC '
                                'SILT_FRAC '
                                'CLAY_FRAC '
                                'OM '
                                'PH_KCL ')
                                )

        profile_df['PROF_NUM'] = profile_df['PROF_NUM'].str[2:]
        profile_df['HOR_THICK'] = profile_df['DEPTHEND'] - profile_df['DEPTHSTART']
        profile_df['NUM_NUMER_COMP'] = \
                profile_df.apply(lambda row: self.num_numeric(row), axis=1)
        profile_df['SAND_FRAC'] = int(profile_df['TEXT_2000']) + int(profile_df['TEXT_200'])
        profile_df['SILT_FRAC'] = int(profile_df['TEXT_20']) + int(profile_df['TEXT_50'])
        profile_df['CLAY_FRAC'] = int(profile_df['TEXT_2'])

        spade_soil_props = spade_data(int(profile_df.PROF_NUM),
                                      int(profile_df.HOR_NUM),
                                      int(profile_df.HOR_THICK),
                                      int(profile_df.NUM_NUMER_COMP),
                                      int(profile_df.SAND_FRAC),
                                      int(profile_df.SILT_FRAC),
                                      int(profile_df.CLAY_FRAC),
                                      float(profile_df.OM),
                                      float(profile_df.PH),
                                      )

        return spade_soil_props

    def param_gen(self, param_array):
        '''
        create a generator for each parameter
        '''
        for line in param_array:
            for item in line:
                yield item

    def rast_extract(self, file_dir, file_name, dtype=int, skip=5):
        '''
        Import the raster file which contains the profile numbers that can
        be linked to the profile names. Store in an efficient numpy array
        '''
        rast_file = file_dir + '\\' + file_name + '.asc'
        with open(rast_file, 'r') as rast_data:
            #skip first five header lines:
            for i in range(skip):
                next(rast_data)
            data_iter = csv.reader(rast_data,
                                   delimiter = ' ')
            #The first character in each line is a space meaning we get an
            #empty string. We don't want this hence data[1:]
            self.data = [data[1:] for data in data_iter]
        rast_array = np.asarray(self.data, dtype=dtype)
        return rast_array


    def hor_perc_calc(self, start, end, row):
        '''
        Calculate the percent of the Toth soil column in the spade columns
        start: start of the Toth soil column
        end:   end of the Toth soil column
        row:   the row of the df
        '''
        if start < row.DEPTHEND:
            if end < row.DEPTHSTART:
                return 0
            if end < row.DEPTHEND:
                if start > row.DEPTHSTART:
                    return 100
                return float(end - row.DEPTHSTART)/(end-start)*100
            return float(row.DEPTHEND - start)/(end-start)*100
        return 0

    def rearrange_dataframe(self):
        '''
        Rearrange the data frame ready for export
        '''
        #first replace weird values with NaN
        self.sol_df.PROF_NUM = self.sol_df.PROF_NUM.apply(float)
        self.sol_df[self.sol_df<-10] = np.NaN

        #fill the NaN values with the data from the previous horizon
        self.sol_df = self.sol_df.fillna(method='pad')

        #Convert the % sand, silt, clay and OM to a fraction
        self.sol_df['SAND_FRAC'] = self.sol_df['SAND_FRAC']/100.0
        self.sol_df['SILT_FRAC'] = self.sol_df['SILT_FRAC']/100.0
        self.sol_df['CLAY_FRAC'] = self.sol_df['CLAY_FRAC']/100.0
        self.sol_df['OM'] = self.sol_df['OM']/100.0

        #Convert Ksat to m/d
        self.sol_df['SAT_HYD_CON'] = self.sol_df['SAT_HYD_CON']/100
        self.sol_df['PHY_SAT_HYD_CON'] = self.sol_df['PHY_SAT_HYD_CON']/100

        #Make sure the following are ints:
        if self.mode == 'HYPRES':
            self.sol_df['PROF_NUM'] = self.sol_df['PROF_NUM'].apply(int)
        else:
            self.sol_df['PROF_NUM'] = self.sol_df.index / 7 + 1
            self.sol_df['PROF_NUM'] = self.sol_df['PROF_NUM'].apply(int)
        self.sol_df['HOR_NUM'] = self.sol_df['HOR_NUM'].apply(int)
        #This is in m and shouldn't be an int? Could be a difference between
        #HYPRES and TOTH:
        #self.sol_df['HOR_THICK'] = self.sol_df['HOR_THICK'].apply(int)/100
        self.sol_df['NUM_NUMER_COMP'] = self.sol_df['NUM_NUMER_COMP'].apply(int)

        #Columns to export in order:
        columns = ['PROF_NUM',
                   'HOR_NUM',
                   'HOR_THICK',
                   'NUM_NUMER_COMP',
                   'SAND_FRAC',
                   'SILT_FRAC',
                   'CLAY_FRAC',
                   'OM',
                   'PH_KCL',
                   'SAT_SOL_WAT_CONT',
                   'RES_WAT_CONT',
                   'PARAM_ALPHA_DRY',
                   'PARAM_ALPHA_WET',
                   'PARAM_N',
                   'SAT_HYD_CON',
                   'PHY_SAT_HYD_CON',
                   'PARAM_L',
                   'DISP_LEN',
                   'SESQ_OX_CONT',
                   ]
        return self.sol_df.loc[:,columns]

    def sol_export(self):
        '''
        Export the data as a .sol file suitable for use in GeoPEARL
        '''
        sol_df_to_export = self.rearrange_dataframe()
        output_file = os.path.join(self.output_folder, "Schematisation.sol")
       
        header = ('*-------------------------------------------------------------------------------\n'
                  '* \n'
                  '* SOIL DATABASE\n'
                  '* =================================\n'
                  '*\n'
                  '* File containing the soil database for the Study area.\n'
                  '* The first part of the file contains paramters that are assumed to be spatially\n'
                  '* constant. The second part of the file contains the spatially distributed\n'
                  '* parameters.\n'
                  '*\n'
                  '*-------------------------------------------------------------------------------\n'
                  '\n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Soil evaporation\n'
                  '\n'
                  'Black      OptSolEvp                    Use the Black option for soil evaporation\n'
                  '0.005      PrcMinEvp          (m.d-1)   Minimum rainfall to reset Black model\n'
                  '0.35       CofRedEvp          (cm1/2)   Reduction parameter in Black equation\n'
                  '1.0        FacEvpSol          (-)       "Crop" factor for soil evaporation\n'
                  '\n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Dispersion length and relative diffusion coefficient\n'
                  '* GeoPEARL only supports the Millington Quirk option!\n'
                  '\n'
                  '2.0               ExpDifLiqMilNom (-)         Exponent in nominator of equation [0.1|5]\n'
                  '0.67              ExpDifLiqMilDen (-)         Exponent in denominator of eqn    [0.1|2]\n'
                  '2.0               ExpDifGasMilNom (-)         Exponent in nominator of equation [0.1|5]\n'
                  '0.67              ExpDifGasMilDen (-)         Exponent in denominator of eqn    [0.1|2]\n'
                  '\n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Depth dependence of transformation\n'
                  'table FacZTra (-)\n'
                  '0.00  1.00\n'
                  '0.30  1.00\n'
                  '0.31  0.50\n'
                  '0.60  0.50\n'
                  '0.61  0.30\n'
                  '1.00  0.30\n'
                  '1.01  0.00\n'
                  '50.00 0.00 \n'
                  'end_table\n'
                  '\n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Depth dependence of sorption\n'
                  'table FacZSor (-)\n'
                  '0.00  1.00\n'
                  '0.30  1.00\n'
                  '0.31  1.00\n'
                  '0.60  1.00\n'
                  '0.61  1.00\n'
                  '1.00  1.00\n'
                  '1.01  1.00\n'
                  '50.00 1.00 \n'
                  'end_table\n'
                  '\n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Column 1  :  Soil profile number\n'
                  '* Column 2  :  Soil horizon number\n'
                  '* Column 3  :  Horizon thickness (m)\n'
                  '* Column 4  :  Number of numerical soil compartments\n'
                  '* Column 5  :  Sand fraction (kg.kg-1) as part of mineral soil\n'
                  '* Column 6  :  Silt fraction (kg.kg-1) as part of mineral soil\n'
                  '* Column 7  :  Clay fraction (kg.kg-1) as part of mineral soil\n'
                  '* Column 8  :  Organic matter content (kg.kg-1)\n'
                  '* Column 9  :  pH-KCl\n'
                  '* Column 10 :  Saturated soil water content (m3.m-3)\n'
                  '* Column 11 :  Residual water content (m3.m-3)\n'
                  '* Column 12 :  Parameter alpha (dry) (cm-1)\n'
                  '* Column 13 :  Parameter alpha (wet) (cm-1)\n'
                  '* Column 14 :  Parameter n (-)\n'
                  '* Column 15 :  Saturated hydraulic conductivity (m.d-1)\n'
                  '* Column 16 :  Physical saturated hydraulic conductivity (m.d-1)\n'
                  '* Column 17 :  Parameter L (-)\n'
                  '* Column 18 :  Dispersion length (m)\n'
                  '* Column 19 :  Sesqui-oxide content (mmol.kg-1)\n'
                  '\n'
                  '\n'
                  '\n'
                  'table SoilProfiles\n'
                  )


        #First round the numbers to a suitable number of decimal places
        round_dict = {'HOR_THICK': 2,
                      'SAND_FRAC': 3,
                      'SILT_FRAC': 3,
                      'CLAY_FRAC': 3,
                      'OM': 3,
                      'PH_KCL': 1,
                      'SAT_SOL_WAT_CONT': 2,
                      'RES_WAT_CONT': 2,
                      'PARAM_ALPHA_DRY': 4,
                      'PARAM_ALPHA_WET': 4,
                      'PARAM_N': 3,
                      'SAT_HYD_CON': 3,
                      'PHY_SAT_HYD_CON': 3,
                      'PARAM_L': 3,
                      'DISP_LEN': 2,
                      }

        sol_df_to_export = sol_df_to_export.round(round_dict)
        sol_df_to_export = sol_df_to_export.sort_values(by = ['PROF_NUM', 'HOR_NUM'])
        with open(output_file, 'w') as sol_file:
            sol_file.write(header)
            #Format for a fixed width text file
            fmt = ('%-8s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s'
                   '%-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s')

            #Need to first convert the df to a numpy array of str so we don't
            #convert int to float
            sol_df_to_export = sol_df_to_export.to_numpy(dtype=str)
            np.savetxt(sol_file, sol_df_to_export, fmt=fmt)
            sol_file.write('end_table\n')

    def plo_export(self):
        '''
        Export the data as a .plo file suitable for use in GeoPEARL
        '''
        output_file = os.path.join(self.output_folder, "Schematisation.plo")

        header = ('*-------------------------------------------------------------------------------\n'
                  '* \n'
                  '* PLOT SCHEMATISATION\n'
                  '* =======================================\n'
                  '*\n'
                  '* File containing the plot schematisation for the Netherlands.\n'
                  '* Plot schematisation is based on Kroon et al. (2001)\n'
                  '*-------------------------------------------------------------------------------\n'
                  ' \n'
                  '*-------------------------------------------------------------------------------\n'
                  '* Column 1  :  The plot ID\n'
                  '* Column 2  :  Area (km2)\n'
                  '* Column 3  :  Meteo district, corresponds with meteo stations in geo file\n'
                  '* Column 4  :  The soil profile number, corresponds with Schemitisation.sol\n'
                  '* Column 5  :  The crop type (1 = grass, 2 = maize, 3 = potato, 4 = nature)\n'
                  '* Column 6  :  Correction factor for precipitation (-)\n'
                  '* Column 7  :  Correction temperature (C)\n'
                  '* Column 8  :  Correction factor for evapotranspiration (-)\n'
                  '* Column 9  :  Irrigation switch\n'
                  '* Column 10 :  Maximum ponding depth (m)\n'
                  '* Column 11 :  Air boundary layer thickness (m)\n'
                  '* Column 12 :  Relative vulnerability rank (1 = lowest score; 2981 = highest score)\n'
                  '* Column 13 :  Texture class of soil profile (1 = sand, 2 = clay, 3 = loam, 4 = peat)\n'
                  '\n'
                  '\n'
                  '\n'
                  'table Plots\n'
                  )
#        self.plo_dict_to_export['vuln_rank'] = self.plo_dict_to_export['plot']
        self.plo_dict_to_export['vuln_rank'] = 1


        with open(output_file, 'w') as plo_file:
            plo_file.write(header)
            #Format for a fixed width text file
            fmt = ('%-13s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s'
                   '%-12s %-12s %-12s %-12s')

            #Need to first convert the df to a numpy array of str so we don't
            #convert int to float
            sol_df_to_export = self.plo_dict_to_export.to_numpy(dtype=str)
            np.savetxt(plo_file, sol_df_to_export, fmt=fmt)
            plo_file.write('end_table\n')

    def map_export(self, raster_type):
        '''
        Export the data as a .map file suitable for use in GeoPEARL
        '''
        output_file = os.path.join(self.output_folder, "Schematisation_uc.map")

        if raster_type == 'HYPRES':
            raster_to_plot = self.hypres_plot_raster
        else:
            raster_to_plot = self.toth_plot_raster

        header = ('ncols        266\n'
                  'nrows        247\n'
                  'xllcorner    -866030.583899999969\n'
                  'yllcorner    -141833.819499999983\n'
                  'cellsize     1000.000000000000\n'
                  )
        with open(output_file, 'w', newline='') as map_file:
            map_file.write(header)
            col_num = 1
            for plot in raster_to_plot:
                if col_num % 266 == 0:
                    map_file.write(str(plot) + '\n')
                    col_num = 1
                else:
                    map_file.write(str(plot) + ' ')
                    col_num += 1
