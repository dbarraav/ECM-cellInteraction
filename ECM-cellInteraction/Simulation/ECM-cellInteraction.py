
from cc3d import CompuCellSetup
        

from ECM-cellInteractionSteppables import ECM-cellInteractionSteppable

CompuCellSetup.register_steppable(steppable=ECM-cellInteractionSteppable(frequency=1))


CompuCellSetup.run()
