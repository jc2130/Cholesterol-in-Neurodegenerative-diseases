import os
import sys

cwd = os.getcwd()

root_folder = os.sep+"team-project"
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep+"utils"+os.sep)

from hfpn import HFPN
from visualisation import Analysis
import AD_HFPN

from validation_data import *
from initial_tokens import *

from datetime import datetime

filenames = ['2x27OH',
             '2xchol',
             'age1',
             'ApoE1',
             'diseased',
             'healthy']


filenames = f'ad_{filenames}_{2million}'

# Confgure and import sub-networks
pn = HFPN(time_step=0.001) # time_step in s
AD_HFPN.add_cholesterol_homeostasis(pn)
AD_HFPN.add_tau_pathology(pn)
AD_HFPN.add_abeta_pathology(pn)
AD_HFPN.add_ER_retraction_collapse(pn)
AD_HFPN.add_energy_metabolism(pn)
AD_HFPN.add_calcium_homeostasis(pn)

### Set diseased condition ###

## 2x 27OH (high cholesterol)##
# disease_multiplier_27OH -> 2 (rather than 1)

## 0.5x 27OH (low cholesterol) ##
# disease_multiplier_27OH -> 0.5 (rather than 1)

## 2x chol ##
#pn.set_initial_tokens_for('p_ApoEchol_extra', it_p_ApoEchol_extra*2)

## 0.5x chol ##
#pn.set_initial_tokens_for('p_ApoEchol_extra', it_p_ApoEchol_extra*0.5)

## age1 ##
#pn.set_initial_tokens_for('p_age', 1) 

## ApoE1 ##
#pn.set_initial_tokens_for('p_ApoE', ad_p_ApoE)


# diseased 
# high 270H and high chol

# healthy -> none



# Run network and plot result
start_time = datetime.now()
pn.run_many_times(1, int(2e6))
execution_time = datetime.now() - start_time

print('\n\ntime to execute:', execution_time)

analysis = Analysis(pn)
Analysis.store_to_file(analysis, )