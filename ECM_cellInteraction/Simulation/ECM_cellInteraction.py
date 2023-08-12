
from cc3d import CompuCellSetup
        

from ECM_cellInteractionSteppables import ECM_cellInteractionSteppable

CompuCellSetup.register_steppable(steppable=ECM_cellInteractionSteppable(frequency=1))


CompuCellSetup.run()
