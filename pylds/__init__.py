"""
PyLDs init
"""

__all__ = ['base', 'config', 'vector_fields', 'tools']

__title__ = "pylds"
__author__ = "Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido"
__copyright__ = "Copyright 2020, PyLDs contributors"
__license__ = "MIT"
__version__ = "0.2.1"
__mail__ = 'ba13026@my.bristol.ac.uk'
__maintainer__ = __author__
__status__ = "Development"


# Default parameters for box-escape condition
# for Variable Time Integration


from .base import *
from .vector_fields import *
from .tools import *
