
from cc3d import CompuCellSetup
        

from ECM_cellInteractionSteppables import ECM_cellInteractionSteppable
from CellAlignment import ECM_AlignmentSteppable

CompuCellSetup.register_steppable(steppable=ECM_cellInteractionSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=ECM_AlignmentSteppable(frequency=1))

CompuCellSetup.run()
