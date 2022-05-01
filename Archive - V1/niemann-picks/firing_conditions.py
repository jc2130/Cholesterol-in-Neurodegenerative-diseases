# Firing conditions for all transitions
from initial_tokens import *
from parameters import *
from rate_functions import *


# Choleserol Homeostasis
fc_t_LDLR_endocyto = lambda a : a["p_ApoEchol_extra"] > 1
fc_t_chol_exocyto = lambda a : a['p_chol_PM'] > 1
fc_t_chol_trans_LE_ER = lambda a : a["p_chol_LE"] > 1
fc_t_chol_trans_LE_mito = lambda a : a["p_chol_LE"] > 1 and a["p_PERK_ER"] > 0
fc_t_chol_trans_LE_PM = lambda a : a["p_chol_LE"] > 1
fc_t_chol_trans_PM_ER = lambda a : a["p_chol_PM"] > 1
fc_t_chol_trans_ER_PM = lambda a : a["p_chol_ER"] > 0
fc_t_chol_trans_ER_mito = lambda a : a["p_chol_ER"] > 0 and a["p_PERK_ER"] > 0
fc_t_CYP11A1_metab = lambda a : a["p_chol_mito"] > 0


# Tau Pathology
fc_t_phos_tau = lambda a : a['p_chol_LE'] >=0 and a['p_tau'] >=0
fc_t_dephos_tauP = lambda a : a['p_tauP'] >=0


# ER Retraction & Collapse
fc_t_RTN3_exp = lambda a: True
fc_t_LE_retro = lambda a : a['p_RTN3_axon'] > 1 and a['p_ATP'] > 1 
fc_t_LE_antero = lambda a : a['p_RTN3_PN'] > 1 and a['p_ATP'] > 1 
fc_t_RTN3_aggregation = lambda a : (a["p_RTN3_axon"] + a["p_RTN3_PN"]) > (it_p_RTN3_tot * (1 - (1-dec_t_RTN3_aggregation) * (a['p_Ab'] - it_p_Ab)/(Ab_t_RTN3_aggregation - it_p_Ab))) # proportional reduction in RTN3 aggregation threshold based on Ab level
fc_t_RTN3_auto = lambda a : a['p_RTN3_HMW_cyto'] > 1 and a['p_RTN3_axon'] > 1
fc_t_RTN3_lyso = lambda a : a['p_RTN3_HMW_auto'] > 1
fc_t_RTN3_dys_auto = lambda a : a['p_RTN3_HMW_cyto'] > 1 and a['p_RTN3_axon'] < it_p_RTN3_axon # any reduction in functional tubular ER from starting value can cause dysfunctional autophagosome formation
fc_t_RTN3_dys_lyso = lambda a : a["p_RTN3_HMW_auto"] > 1


# Calcium and Sphingosine Homeostasis
fc_t_RyR_IP3R = lambda a : a['p_Ca_extra'] == 0
fc_t_Ca_cyto_LE = lambda a : a['p_Ca_cyto'] >= 1
fc_t_TRPML1 = lambda a : a['p_Ca_LE'] >= 1
fc_t_Ca_imp = lambda a : a['p_Ca_extra'] == 1
fc_t_NCX_PMCA = lambda a: a['p_on3'] == 1 
fc_t_mNCLX = lambda a : a['p_Ca_mito'] >0 
fc_t_SERCA = lambda a : a['p_ATP'] > 0 
fc_t_mCU = lambda a : a['p_Ca_cyto'] > it_p_Ca_cyto
fc_t_MAM = lambda a : a['p_Ca_ER'] > it_p_Ca_cyto 
fc_t_NaK_ATPase = lambda a: a['p_on3'] == 1 and a['p_ATP'] > 0


# Mitochondria
fc_t_krebs = lambda a : a['p_ADP'] > 0
fc_t_ATP_hydro_mito = lambda a : a['p_ATP'] > 0
fc_t_ETC = lambda a : a['p_reduc_mito'] > 0 and a['p_ADP'] > 0
fc_t_ROS_metab = lambda a : a['p_ROS_mito'] > 0 


# Apoptosis / Survival
fc_t_ERstress_Ca_cyto_induc = lambda a : a['p_ROS_mito'] >=0 and a['p_tauP'] >=0
fc_t_BAX_actv = lambda a : a['p_PERK_ER'] >=0 and a['p_Bcl2_mito'] >= 0
fc_t_Bcl2_actv = lambda a : a['p_PERK_ER'] >=0
fc_t_mito_dysfunc = lambda a : a['p_BAX_mito'] >=0 and a['p_ROS_mito'] >=0
fc_t_TRADD_actv = lambda a : a['p_chol_LE'] >=0
fc_t_Cas3_actv = lambda a : a['p_TRADD_cyto'] >=0 and a['p_cytc_cyto'] >=0
fc_t_Cas3_degr = lambda a : a['p_cas3'] >= 0
fc_t_mTORC1_actv = lambda a : a['p_Ca_LE'] >=0 and a['p_NPC1_LE'] >=0
fc_t_NFkB_actv = lambda a : a['p_TRADD_cyto'] >= 0 or a['p_mTORC1_LE'] >=0
fc_t_NFkB_degr = lambda a : a['p_NFkB_cyto'] > 0