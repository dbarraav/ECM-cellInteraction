from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *
import numpy as np
import os 
import csv

totalSimTime = 1000
magn_ECM = 1.0
force = 10.0

alignment_time = 10

class ECM_AlignmentSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
                
        ECM_A_Vect = np.random.uniform(-np.pi, np.pi, size = len(self.cell_list_by_type(self.ECM)))

        print(ECM_A_Vect)
        ECM_A = self.create_vector_field_cell_level_py("ECM_A")
        
        for val, cell in enumerate(self.cell_list_by_type(self.ECM)):
            ECM_A[cell] = [magn_ECM * np.cos(ECM_A_Vect[val]),
                           magn_ECM * np.sin(ECM_A_Vect[val]), 
                           0]
                           
                           
        for cell in self.cell_list_by_type(self.EPI, self.MES): 
            cell.dict['pos_alignment_time_ago'] = [cell.xCOM, cell.yCOM, cell.zCOM]
        
    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        
        # for each cell, find the average ECM below
        ECM_A_field = self.field.ECM_A

        
        for cell in self.cell_list_by_type(self.MES): 
            
            # get average ECM field below cell
            vector_x = 0
            vector_y = 0
            
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor and neighbor.type == self.ECM:
                    vector_x += ECM_A_field[neighbor][0]
                    vector_y += ECM_A_field[neighbor][1]
                                
            # set direction of active force in the direction of mean ECM direction
            angle = np.arctan2(vector_y, vector_x)
            
            cell.dict['force_angle'] = angle
            
            lambdaVecX = force * np.cos(angle)
            lambdaVecY = force * np.sin(angle)
                    
                    
         
        if mcs%alignment_time == 0: 
            
            for cell in self.cell_list_by_type(self.MES):
                angle_velocity = np.arctan2(cell.yCOM - cell.dict['pos_alignment_time_ago'][1],
                                            cell.xCOM - cell.dict['pos_alignment_time_ago'][0])
          
                # realign the ECM according to the direction of MOTION (!) of the cell (TEST!)
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.ECM:
                        ECM_A_field[neighbor] = [np.cos(angle_velocity), np.sin(angle_velocity),0]
                        
                cell.dict['pos_alignment_time_ago'] = [cell.xCOM, cell.yCOM, cell.zCOM]
                        
     
     
        if mcs >= totalSimTime:
            self.stop_simulation()


        