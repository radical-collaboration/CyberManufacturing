#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 12:32:29 2017

@author: Chaitanya Sampat
"""

'''
This class extracts the data from the liggghts output file at a specified timestep and 
stores the data into dictionaries that can be used for further calculations.
'''
import numpy as np
import os
from pathlib2 import Path

class Controller_DEM_resource_reader(object):
# Class variables / global variables

# Constructor
    def __init__(self, timestep, types, liggghts_output_path):
        self.type = types
        self.raw_collision_data = '';
        self.raw_impact_data = '';
        self.timestep = timestep;
        self.collision_data_acc_types = {};
        self.number_of_impacts = 0.0
        self.number_of_particles = 0.0
        self.liggghts_output_dir = liggghts_output_path
        self.impact_data_acc_types = {};
        self.fields_collision = {'type':[],'velocity':[],'collision':[],'diameter':[]}
        self.fields_impact =  {'velocity':[],'ids':[],'force':[],'contactarea':[],'delta':[]}
        for i in range(1,(types+1)):
            self.collision_data_acc_types.update({str(i):[]})
        self.a = self.initial_all_functions()


    def initial_all_functions(self):
        self.raw_collision_data = self.liggghts_collision_raw_data()
        self.raw_impact_data = self.liggghts_impact_raw_data()
        self.collision_data_acc_types = self.liggghts_collision_data_store()
        self.number_of_impacts = self.liggghts_impact_data_store()


# function to read the certain liggghts present in the folder and return liggghts collision data
    def liggghts_collision_raw_data(self):
        # check if both impact and the collision files exist
        liggghts_output_files_path = self.liggghts_output_dir
        timestep_collision_file = "collision%d.atom"% self.timestep
        timestep_impact_file = 'impact%d.atom'% self.timestep
        path_collision_file = liggghts_output_files_path + str(timestep_collision_file)
        path_impact_file = liggghts_output_files_path + str(timestep_impact_file)
        temp1 = Path(path_collision_file)
        temp2 = Path(path_impact_file)
        if (temp1.is_file() and temp2.is_file()):
            print "Both impact and collision files of %d timestep exists"% self.timestep
        else:
            print "collision file exists at %d : %s"%(self.timestep, str(temp1.is_file()))
            print "impact file exists at %d : %s"%(self.timestep,str(temp2.is_file()))
        collision_file_open = open(path_collision_file, "r")
        #impact_file_open = open(path_impact_file, "r")
        #print file_open.read()
        self.raw_collision_data = collision_file_open.readlines()
        #self.raw_impact_data = impact_file_open.readlines()
        collision_file_open.close()
        #impact_file_open.close()
        return self.raw_collision_data

# function to read the certain liggghts present in the folder and return liggghts collision data
    def liggghts_impact_raw_data(self):
        liggghts_output_files_path = self.liggghts_output_dir
        timestep_impact_file = 'impact%d.atom'% self.timestep
        path_impact_file = liggghts_output_files_path + str(timestep_impact_file)
        impact_file_open = open(path_impact_file, "r")
        self.raw_impact_data = impact_file_open.readlines()
        impact_file_open.close()
        return self.raw_impact_data


# function to store the collision data into a nested dictionary
    def liggghts_collision_data_store(self):
        raw_col_data = self.raw_collision_data
        self.number_of_particles = self.raw_collision_data[3].split(' ')[0]
        temp_col = []
        temp_fields = self.fields_collision
        data_temp = self.collision_data_acc_types
#store data in dictionary according to type
        for x in xrange(0,int(self.number_of_particles)):
            temp_fields['collision'] = []
            temp_fields['velocity'] = []
            temp_fields['diameter'] = []
            temp_fields['type'] = []
            temp_type = raw_col_data[x+9].split(' ')[1]
            temp_fields['type'] = temp_type
            temp_vel = [float(raw_col_data[x+9].split(' ')[8]),float(raw_col_data[x+9].split(' ')[9]),float(raw_col_data[x+9].split(' ')[10])]
            for i in xrange(1,self.type+1):
                temp_col.append(int(raw_col_data[x+9].split(' ')[(i+13)]))
            temp_diam = float(raw_col_data[x+9].split(' ')[31])
            temp_fields['velocity'].append(temp_vel)
            temp_fields['collision'].append(temp_col)
            temp_fields['diameter'].append(temp_diam)
            data_temp[temp_type].append(temp_fields)
            self.collision_data_acc_types[temp_fields['type']].pop()
            self.collision_data_acc_types[temp_fields['type']].append(temp_fields.copy())
            temp_diam = []
            temp_type = []
            temp_col = []
            temp_fields.clear()
<<<<<<< HEAD

        print(len(self.collision_data_acc_types['6']))




        '''
        for i in enumerate(file_open):
            if i == 4:
                number_of_particles = float(linecache.getline(path_collision_file, i))
                number_of_particles = file_open.readlines()[i - 1]

        print(number_of_particles)
        '''


# ------------------------------------------------------------------------------

# abcd = Controller_DEM_resource_reader(200000,16)
# ad = abcd.liggghts_raw_collision_data()
# ad = abcd.ligghts_data_store()
#ad.read()
=======
        return self.collision_data_acc_types

#function to store impact data into a dictionary
    def liggghts_impact_data_store(self):
        raw_imp_data = self.raw_impact_data
        self.number_of_impacts = self.raw_impact_data[3].split(' ')[0]
        temp_imp_fields = self.fields_impact
        for x in xrange(0,int(self.number_of_impacts)):
            temp_imp_fields['velocity'] = []
            temp_imp_fields['ids'] = []
            temp_imp_fields['force'] = []
            temp_imp_fields['contactarea'] = []
            temp_imp_fields['delta'] = []
            temp_imp_vel = [float(raw_imp_data[x+9].split(' ')[0]),float(raw_imp_data[x+9].split(' ')[1]),float(raw_imp_data[x+9].split(' ')[2]),float(raw_imp_data[x+9].split(' ')[3]),float(raw_imp_data[x+9].split(' ')[4]),float(raw_imp_data[x+9].split(' ')[5])]
            temp_imp_ids = [float(raw_imp_data[x+9].split(' ')[6]),float(raw_imp_data[x+9].split(' ')[7]),float(raw_imp_data[x+9].split(' ')[8])]
            temp_imp_force = [float(raw_imp_data[x+9].split(' ')[9]),float(raw_imp_data[x+9].split(' ')[10]),float(raw_imp_data[x+9].split(' ')[11])]
            temp_imp_contactarea = float(raw_imp_data[x+9].split(' ')[12])
            temp_imp_delta = float(raw_imp_data[x+9].split(' ')[13])
            temp_imp_fields['velocity'].append(temp_imp_vel)
            temp_imp_fields['ids'].append(temp_imp_ids)
            temp_imp_fields['force'].append(temp_imp_force)
            temp_imp_fields['contactarea'].append(temp_imp_contactarea)
            temp_imp_fields['delta'].append(temp_imp_delta)
            self.impact_data_acc_types.update({str(x):[]})
            #self.impact_data_acc_types.pop()
            self.impact_data_acc_types[str(x)].append(temp_imp_fields.copy())
            temp_imp_vel = []
            temp_imp_ids = []
            temp_imp_force = []
            temp_imp_contactarea = 0.0
            temp_imp_delta = 0.0
            temp_imp_fields.clear()
        return self.number_of_impacts

# ------------------------------------------------------------------------------
>>>>>>> devel
