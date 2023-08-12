
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

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        # ECM_D = self.field.ECM_D
        ECM_D = self.get_field_secretor("ECM_D")
        for cell in self.cell_list:
            
            # value = ECM_D[cell.xCOM, cell.yCOM, 0]
            cell.dict["ECM_DVal"] = ECM_D.amountSeenByCell(cell)
            
            # (self.EGF_GrowthScalar_STEM) * (cell.dict['EGF']**4/(self.EGF_STEM_HalfMaxValue**4 + cell.dict['EGF']**4)))
            
            H = cell.dict["ECM_DVal"]**n/((thres*cell.volume)**n + cell.dict["ECM_DVal"]**n)
            
            if H > thres:
                cell.type = self.MES
        #

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return


        