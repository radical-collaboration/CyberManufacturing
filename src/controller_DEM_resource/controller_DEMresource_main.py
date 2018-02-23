#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 09:45:37 2017

@author: Chaitanya Sampat
"""

'''
This class is the main execution file for the DEM controller. It needs initial time and number of type of
particles as inputs. This calls the relevent methods from the interpretor and reader class. It starts 
execution at the inital timestep and then checks for the existence of the file at a predefined timestep
difference (defined as dump_difference). Once the files are found then it calls and the other 2 classes
and checks for the output status from the interpretor class
 ========================================================================================================
 OUTPUT STATUS
 
 0 = continue the execution of the DEM, the % change is less than the defined limit
 1 = kill DEM and execute PBM since the change is  more than 10%
 2 = kill DEM and report that the system has reached a steady state and stop all executions.
 ========================================================================================================
''' 
import numpy as np
import os
import time
import controller_DEMresource_data_reader as DEMreader
import controller_DEMresource_data_interpretor as DEMinter
import json
import sys

class controller_DEMresource_main(object):
    def __init__(self, init_timestep, types, liggghts_output_path, pbm_timestep, DEMInitTimeStep, stepRestart, min_dia, max_dia, total_flow_rate, solid_density):
        # initial timestep and type of partilces are taken as input
        self.init_timestep = init_timestep 
        self.type = types
        self.dump_difference = 50000 # this is the predefined interval after which LIGGGHTS ouputs a dump file. Has to be changed if changed in the LIGGGHTS input file
        self.liggghts_output_dir = liggghts_output_path
        self.pbm_timestep = pbm_timestep
        self.DEMInitTimeStep = DEMInitTimeStep
        self.stepRestart = stepRestart
        self.min_dia = min_dia
        self.max_dia = max_dia
        self.total_flow_rate = total_flow_rate
        self.solid_density = solid_density

# ------------------------------------------------------------------------------------------------
    def main(self):
        # This method executes the necessary commands to run the DEM controller and waits for the command from the interpertor to output certain files that the executioner will act upon.
        timestep = int(self.init_timestep) + int(self.dump_difference)
        flag = 0
        dump_data = {}
        timestep_collision_file = "/collision%d.atom"%int(timestep)
        timestep_impact_file = '/impact%d.atom'%int(timestep)
        collision_file = self.liggghts_output_dir + str(timestep_collision_file)
        impact_file = self.liggghts_output_dir + str(timestep_impact_file)
        # the execution waits for a seconds everytime it enters the loop till the file do not exist
        while not(os.path.exists(collision_file) and os.path.exists(impact_file)):
            time.sleep(5)
            print("Waiting for file to be printed")
        # opening the different files so as to write the calculated average / total data ovetime
        dump_avg_vel = open("velocity_average_overtime_starting_%d.txt"%self.init_timestep, "w") 
        dump_collisions = open("collisions_overtime_starting_%d.txt"%self.init_timestep, "w")
        dump_impacts = open("impacts_overtime_starting_%d.txt"%self.init_timestep, "w")
        obj_inter = DEMinter.controller_DEM_resource_interpretor(timestep, self.type, timestep, self.liggghts_output_dir)
        # Keeps checking for the existence of the file till one of the criteria for a killing the DEM are not met
        print 'start searching'
        count = 0
        while (flag == 0 and count < 40):
            # defining the files and path of the files that it needs to search for.
            #liggghts_output_files_path = os.getcwd() + '/liggghts_output_files/'

            timestep_collision_file = "/collision%d.atom"%timestep
            timestep_impact_file = '/impact%d.atom'%timestep
            collision_file = self.liggghts_output_dir + str(timestep_collision_file)
            impact_file = self.liggghts_output_dir + str(timestep_impact_file)
            print(impact_file)
            # the execution waits for a seconds everytime it enters the loop till the file do not exist
            while not(os.path.exists(collision_file) and os.path.exists(impact_file)):
                time.sleep(6)
                count = count + 1
                print("Waiting for file to be printed")
            #once the files are found the it performs the comparison to the initial data files using the interpretor class
            if (os.path.isfile(collision_file) and os.path.isfile(impact_file)):
                flag = obj_inter.liggghts_data_comparison(timestep)
                a = obj_inter.avg_all_data(timestep)
                dump_avg_vel.write("%d %s \n"%(timestep, str(obj_inter.avg_vel_array)))
                dump_collisions.write("%d %d \n"%(timestep, sum(sum(obj_inter.collision_matrix))))
                dump_impacts.write("%d %d \n"%(timestep, int(obj_inter.avg_impacts)))
                timestep = timestep + self.dump_difference
            else :
                continue
        
        dump_avg_vel.close()
        dump_collisions.close()
        dump_impacts.close()

        if (flag == 0):
            statusfile['status'] = []
            statusfile['status'].append(str(flag))

            with open('DEM_status.json' , 'w') as demsf:
                json.dump(status, demsf)

            if stepRestart:
                final_timestep = timestep + 200000
                statusfile['DEMInitTimeStep'].append(self.DEMInitTimeStep)
            else:
                final_timestep = timestep + 500000
                statusfile['DEMInitTimeStep'].append(self.init_timestep)
            # defining parameters for the liggghts restart file
            liggghts_restart_file = liggghts_restart.liggghts_input_creator(self.init_timestep, final_timestep, \
                                        self.min_dia, self.max_dia, self.type, self.total_flow_rate, self.solid_density)
            liggghts_restart_file.main_writer()

        elif (flag == 1):
            # kill the DEM and start the PBM and also print the input file for the PBM executable
            print("Time to change to PBM and kill DEM")
            tot_part_each_type = np.zeros(16)
            avg_vel_array = np.zeros_like(tot_part_each_type)
            tot_part_each_type = obj_inter.tot_part_each_type
            avg_vel_array = obj_inter.avg_vel_array
            tot_part_each_type.astype(int)
            avg_vel_array.astype(float)
            mixing_time = 15.0
            if (self.pbm_timestep > 0):
                mixing_time = 0.0
            timestep = timestep - self.dump_difference
            # Dump the data relevant for the PBM execution as well as the for the DEM restart (if required) into a json file
            dump_data['types'] = []
            dump_data['types'].append(self.type)
            dump_data['total number of particles'] = []
            dump_data['total number of particles'].append(int(obj_inter.num_of_particles))
            dump_data['number of particles of each type'] = []
            dump_data['number of particles of each type'].append(list(tot_part_each_type))
            dump_data['average velocity of each type'] = []
            dump_data['average velocity of each type'].append(list(avg_vel_array))
            dump_data['DEM_time_step'] = []
            dump_data['DEM_time_step'].append(int(timestep))
            dump_data['Aggregation Kernel Constant'] = []
            dump_data['Aggregation Kernel Constant'].append(1e-9)
            dump_data['Breakage Kernel Constant'] = []
            dump_data['Breakage Kernel Constant'].append(1e-7)
            dump_data['mixing_times'] = []
            # change to actual mixing time before running on hpc
            dump_data['mixing_times'].append(15)
            dump_data['Liquid Addition Time'] = []
            dump_data['Liquid Addition Time'].append(45)
            dump_data['PBM_init_time_step'] = []
            dump_data['PBM_init_time_step'].append(self.pbm_timestep)
            #dump_data['collision_matrix'] = []
            #dump_data['collision_matrix'].append(list(list(obj_inter.collision_matrix)))
            #dump_data.update({'collision_matrix': obj_inter.collision_matrix})
            # print(dump_data)
            # status = {'status':str(flag)}
            dump_data['status'] = []
            dump_data['status'].append(str(flag))
            with open('DEM_status.json' , 'w') as demsf:
                json.dump(dump_data, demsf)
            # with open('PBM_input.json' , 'w') as outfile:
            #     json.dump(dump_data, outfile)
            # # also printing a text file, whichever is easier for yuktesh to read, we can remove the json or this once decided
            with open('PBM_input.in', 'w') as ipt:
                # ipt.write("types %d"%int(obj_inter.num_of_particles))
                # ipt.writelines("total number of particles %f\n"%item for item in tot_part_each_type)
                # ipt.writelines("average velocity of particles  %f\n"%item for item in avg_vel_array)
                ipt.write("#PBM input parameters")
                ipt.write("\nAggregationKernelConstant %s"%(str(1.0e-4)))
                ipt.write("\nBreakageKernelConstant %s"%(str(5.0e-20)))
                ipt.write("\nConsolidationConstant %s"%(str(1.0e-20)))
                ipt.write("\nInitialPorosity %f"%(0.4))
                ipt.write("\nMinimumPorosity %f"%(0.2))
                ipt.write("\nGranuleSaturationFactor %f"%(0.5))
                ipt.write("\nCriticalStokesDeformationNumber %f"%(5.9e-3))
                ipt.write("\nBinderViscosity %f"%(0.05))
                ipt.write("\nCoefficientOfRestitution %f"%(0.2))
                ipt.write("\nLiquidThickness %f"%(1.0))
                ipt.write("\nSurfaceAsperities %f"%(5.0e-6))
                ipt.write("\nGranulatorLength %f"%(0.38))
                ipt.write("\nNumberOfCompartmets %d"%(8))
                ipt.write("\nParticleResidenceTime %f"%(5.0))
                ipt.write("\nDiammeterOfImpeller %f"%(0.114))
                ipt.write("\nSpeedOfImpeller %f"%(2000.0))
                ipt.write("\nPreMixingTime %f"%mixing_time)
                ipt.write("\nLiquidAdditionTime %f"%(45.0))
                ipt.write("\nPostMixingTime %f"%(0.0))
                ipt.write("\nThroughput %f"%(15.0))
                ipt.write("\nSolidDensity %f"%(476.0))
                ipt.write("\nLiquidToSolidRatio %f"%(0.35))
                ipt.write("\nLiquidDensity %f"%(1000.0))
                ipt.write("\nNumberOfFirstSolidBins %d"%(16))
                ipt.write("\nFirstSolidVolumeCoefficients %s"%(str(5.0e-16)))
                ipt.write("\nFirstSolidVolumeBase %f"%(3.0))
                ipt.write("\nNumberOfSecondSolidBins %d"%(16))
                ipt.write("\nSecondSolidVolumeCoefficient %s"%(str(5.0e-16)))
                ipt.write("\nSecondSolidVolumeBase %f"%(3.0))
                ipt.write("\nCollisionEfficiencyConstant %f"%(0.01))
                ipt.write("\nDEMTimeStep %s"%(str(5.0e-7)))
                ipt.write("\nNumberOfDEMBins %d"%(16))

        elif (flag == 2):
            # kill the DEM since there has been no change in the number of collisions / impacts / velocity for 2 seconds.
            status = {'status':str(flag)}
            status['time step'] = []
            status['time step'].append(int(timestep))
            with open('DEM_status.json' , 'w') as demsf:
                json.dump(status, demsf)
            print("The system is at steady state")


# ------------------------------------------------------------------------------------------------

abcd = controller_DEMresource_main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6], float(sys.argv[7]), \
                                        float(sys.argv[8]), float(sys.argv[9], float(sys.argv[10]))
abcd.main()
