
from cc3d.core.PySteppables import *
n  = 2
thres = 0.5

class ECM_cellInteractionSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        
        
        for x in range(0, self.dim.x, 4):
            for y in range(0, self.dim.y, 4): 
                self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.ECM)
                
        for x in range(0, self.dim.x - 7, 7):
            for y in range(0, self.dim.y - 7,7):
                if (x-self.dim.x/2)**2 + (y-self.dim.y/2)**2 < 68**2: 
                    self.cell_field[x:x+7, y:y+7, 1] = self.new_cell(self.EPI)
                    
    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        # ECM_D = self.field.ECM_D
        ECM_D = self.get_field_secretor("ECM_D")
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
            
            
            # self.adhesionFlexPlugin.setMediumAdhesionMoleculeDensityByIndex(0, 11.2)
                
        

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return


        