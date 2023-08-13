from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *
import numpy as np
n  = 2
thres = 0.5

# ECM_A_Vals = np.random.

class ECM_cellInteractionSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        self.create_scalar_field_cell_level_py("molCM")
        
        for x in range(0, self.dim.x, 4):
            for y in range(0, self.dim.y, 4): 
                self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.ECM)
                
        for x in range(0, self.dim.x - 7, 7):
            for y in range(0, self.dim.y - 7,7):
                if (x-self.dim.x/2)**2 + (y-self.dim.y/2)**2 < 68**2: 
                    self.cell_field[x:x+7, y:y+7, 1] = self.new_cell(self.EPI)
                    
        # ECM_A_Vals = np.random.rand(len(self.cell_list_by_type(self.ECM)))
        
        # for val, cell in enumerate(self.cell_list_by_type(self.ECM)):
            # cell.dict["ECM_A_Val"] = ECM_A_Vals[val]
            # # Make sure Secretion plugin is loaded
            # # make sure this field is defined in one of the PDE solvers
            # # you may reuse secretor for many cells. Simply define it outside the loop
            # secretor = self.get_field_secretor("ECM_A")
            # secretor.secreteOutsideCellAtBoundary(cell, cell.dict["ECM_A_Val"])           
        
        self.create_vector_field_cell_level_py("ECM_A")
            
        
            
    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        CMField = self.field.molCM
        # CMField[:,:,:] = 0
        # ECM_D = self.field.ECM_D
        ECM_D = self.get_field_secretor("ECM_D")
        ECM_A = self.get_field_secretor("ECM_A")
        for cell in self.cell_list_by_type(self.EPI):
            cellECMContact = 0 
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if not neighbor:
                    cellECMContact = common_surface_area

            # value = ECM_D[cell.xCOM, cell.yCOM, 0]
            cell.dict["ECM_DVal"] = ECM_D.amountSeenByCell(cell)
            
            # (self.EGF_GrowthScalar_STEM) * (cell.dict['EGF']**4/(self.EGF_STEM_HalfMaxValue**4 + cell.dict['EGF']**4)))
            
            H = cell.dict["ECM_DVal"]**n/((thres*cell.volume)**n + cell.dict["ECM_DVal"]**n)
            
            if (H > thres) and cellECMContact > 0:
                cell.type = self.MES
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'CC',5.5)
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'CM',3)
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'MM',0)
        
        for cell in self.cell_list_by_type(self.MES):
            H = cell.dict["ECM_DVal"]**n/((thres*cell.volume)**n + cell.dict["ECM_DVal"]**n)
            if (H > thres):
                HH = 10*H
                # print(HH)
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'CC',5.5)
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'CM', HH)
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell,'MM',0)
                
            CMField[cell] = self.adhesionFlexPlugin.getAdhesionMoleculeDensity(cell, "CM")
            
            maxSurf = 0
            migDir = np.array([0, 0 , 0])
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                
                if neighbor:
                    if neighbor.type == self.ECM:
                        # alignDetect = ECM_A.amountSeenByCell(neighbor)
                        alignDetect = neighbor.dict["ECM_A_Val"]
                        if alignDetect > maxSurf:
                            maxSurf = alignDetect
                            neighborCOM = np.array([neighbor.xCOM, neighbor.yCOM])
                        # print(neighbor)
                        
            # print(neighborCOM)
            shift_vector = neighborCOM - np.array([cell.xCOM, cell.yCOM])
            shift_vector = shift_vector/np.linalg.norm(shift_vector)
            # self.move_cell(cell, shift_vector)
            cell.lambdaVecX = - 200 * shift_vector[0] # force component pointing along X axis - towards positive X's
            cell.lambdaVecY = - 200 * shift_vector[1] # force component pointing along Y axis - towards negative Y's
            # cell.lambdaVecZ = 0.0  force component pointing along Z axis
                        
                        
                
            
            
                    
            
            
            alignDetect = ECM_A.amountSeenByCell(cell)
            # print('this is the cell id {} and this is fiberAlignment {}'.format(cell.id, alignDetect))
            
            # self.adhesionFlexPlugin.setMediumAdhesionMoleculeDensityByIndex(0, 11.2)
                
        

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return


        