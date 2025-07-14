import numpy as np
import sunpy.map
import sys
from scc_funs import scc_make_array, scc_zelensky_array
from cor_prep import cor_prep
from hi_prep import hi_prep
from astropy.io import fits

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)
np.seterr(divide='ignore')


def wispr_prep():
    print ('hi')
    return 6, 5