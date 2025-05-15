"""
Math functions and further calculations.
"""

from pencil.util import pc_print
from .streamlines import *

####### working with particles and grid data
from .part_to_grid import *  # bin particle quantities to a grid
from .fill_gaps_in_grid import *
from .accuracy import *
from .draglift import *
from .gravity import grav_profile
from .tensors import *
from .Reynolds import *
from .shocktube import sod
from .Gaussian_averages import kernel_smooth, gauss_3Dsmooth

try:
    from .aver2h5 import *
except:
    pc_print("Warning: Could not import calc.aver2h5. Try:")
    pc_print("'pip3 install h5py' (Python 3) or 'pip install h5py' (Python 2).")
