# Initial tokens of all places
from parameters import *

# Choleserol Homeostasis
it_p_ApoEchol_extra = 1.60*10**6
it_p_chol_LE = 2.65*10**9            
it_p_chol_mito = 3.79*10**9      
it_p_chol_ER = 4.37*10**8
it_p_chol_PM = 8.45*10**10
it_p_preg = 0
it_p_NPC1_LE = 1 #660000 for diseased state
it_p_SINK = 1


# Tau Pathology
it_p_GSK3b_inact = 2684994
it_p_GSK3b_act = 2815005
it_p_tau = 3.72*10**7
it_p_tauP = 6.36*10**6
it_p_PP_2A = 1.33e7


# ER Retraction & Collapse
it_p_RTN3_axon =   5344426
it_p_RTN3_PN =   6655574
it_p_RTN3_tot = it_p_RTN3_axon + it_p_RTN3_PN 
it_p_RTN3_HMW_cyto = 0
it_p_RTN3_HMW_auto = 0
it_p_RTN3_HMW_lyso = 0
it_p_RTN3_HMW_dys1 = 0
it_p_RTN3_HMW_dys2 = 0
it_p_Ab = 34048  


# Calcium and Sphingosine Homeostasis
it_p_Ca_LE = 2.44*10**5
it_p_Ca_extra = 1
it_p_Ca_cyto = 7.2*10**5
it_p_Ca_ER = 1.8*10**9
it_p_Ca_mito = 2968656.262
it_p_LRRK2_mut = 0 # PD specific token that is used in r_t_mNCLX 

# Mitochondria
it_p_ATP = 5.42*10**9
it_p_ADP = 5.42*10**6
it_p_reduc_mito = 2.71*10**9
it_p_ROS_mito = 13557.04
it_p_H2O_mito = 0

# Apoptosis / Survival
it_p_PERK_ER = 1*10**6  #Do NOT change this anymore for disease state, only NPC1_LE to 660000!
it_p_BAX_mito = 0.972
it_p_Bcl2_mito = 30000
it_p_cytc_cyto = 13015
it_p_TRADD_cyto = 298540
it_p_mTORC1_LE = 5400
it_p_cas3 = 40000 
it_p_NFkB_cyto = 100