import numpy as np
import os

# class to read the ligghts data file, process data and create PBM input file

class Controller_DEM_resource_reader(object):
# Class variables / global variables

# Constructor
    def __init__(self, timestep, types):
        self.type = types
        self.raw_collision_data = '';
        self.timestep = timestep;
        self.collision_data_acc_types = {};
        self.fields={'type':[],'velocity':[],'collision':[],'diameter':[]}
        for i in range(1,(types+1)):
            self.collision_data_acc_types.update({str(i):[]})


# function to read the certain liggghts present in the folder and return data
    def liggghts_raw_collision_data(self):
        liggghts_output_files_path = os.getcwd() + '/liggghts_output_files/'
        timestep_collision_file = "collision%d.atom"% self.timestep
        path_collision_file = liggghts_output_files_path + str(timestep_collision_file)
        collision_file_open = open(path_collision_file, "r")
        #print file_open.read()
        self.raw_collision_data = collision_file_open.readlines()
        collision_file_open.close()
        return self.raw_collision_data

# function to store the data into a nested dictionary
    def ligghts_data_store(self):
        raw_col_data = self.raw_collision_data
        number_of_particles = self.raw_collision_data[3].split(' ')[0]
        temp_col = []
        temp_fields = self.fields
#print number_of_particles
        for x in xrange(0,int(number_of_particles)):
            data_temp = self.collision_data_acc_types
            temp_fields['collision'] = []
            temp_fields['velocity'] = []
            temp_fields['diameter'] = []
            temp_fields['type'] = []
            temp_type = raw_col_data[x+9].split(' ')[1]
            temp_fields['type'] = temp_type
            temp_vel = [float(raw_col_data[x+9].split(' ')[8]),float(raw_col_data[x+9].split(' ')[9]),float(raw_col_data[x+9].split(' ')[10])]
            for i in xrange(1,17):
                temp_col.append(int(raw_col_data[x+9].split(' ')[(i+13)]))
            temp_diam = float(raw_col_data[x+9].split(' ')[31])
            #print (temp_fields)
            temp_fields['velocity'].append(temp_vel)
            temp_fields['collision'].append(temp_col)
            temp_fields['diameter'].append(temp_diam)
            data_temp[temp_type].append(temp_fields)
            print(temp_fields['type'])
            #self.collision_data_acc_types[temp_fields['type']].append(temp_fields)
            temp_diam = []
            temp_type = []
            temp_col = []
            #temp_fields.clear()
        print(self.collision_data_acc_types['5'])




        '''
        for i in enumerate(file_open):
            if i == 4:
                number_of_particles = float(linecache.getline(path_collision_file, i))
                number_of_particles = file_open.readlines()[i - 1]

        print(number_of_particles)
        '''


# ------------------------------------------------------------------------------

abcd = Controller_DEM_resource_reader(200000,16)
ad = abcd.liggghts_raw_collision_data()
ad = abcd.ligghts_data_store()
#ad.read()
