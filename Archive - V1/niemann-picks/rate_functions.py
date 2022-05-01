# Rate functions for all transitions
from parameters import *
from initial_tokens import *

# Choleserol Homeostasis
r_t_LDLR_endocyto = lambda a : it_p_chol_ER**(0.5) / a["p_chol_ER"]**(0.5) * vmax_t_LDLR_endocyto * min([1, max([0.3,(m_t_LDLR_endocyto * (a["p_chol_ER"]/it_p_chol_ER) + n_t_LDLR_endocyto)])]) * a["p_ApoEchol_extra"]/(Km_t_LDLR_endocyto + a["p_ApoEchol_extra"])
r_t_chol_trans_LE_ER = lambda a : k_t_chol_trans_LE_ER * a["p_chol_LE"] / a["p_NPC1_LE"]
r_t_chol_trans_LE_mito = lambda a : k_t_chol_trans_LE_mito * a["p_chol_LE"] * a["p_PERK_ER"]**(1/2)
r_t_chol_trans_LE_PM = lambda a : k_t_chol_trans_LE_PM * a["p_chol_LE"]
r_t_chol_trans_ER_PM = lambda a : k_t_chol_trans_ER_PM * a["p_chol_ER"]
r_t_chol_trans_PM_ER = lambda a : k_t_chol_trans_PM_ER * a["p_chol_PM"]
r_t_chol_trans_ER_mito = lambda a : k_t_chol_trans_ER_mito * a["p_chol_ER"] * a["p_PERK_ER"]**(1/2)
r_t_CYP11A1_metab = lambda a : ((vmax_t_CYP11A1_metab * a["p_chol_mito"])/(Km_t_CYP11A1_metab + a["p_chol_mito"]) + (vmax_t_CYP27A1_metab* a["p_chol_mito"])/(Km_t_CYP27A1_metab + a["p_chol_mito"]))*a['p_chol_mito']/it_p_chol_mito


# Tau Pathology
r_t_actv_GSK3b = lambda a : tau_multiplier * k_t_actv_GSK3b * a['p_GSK3b_inact'] * (dis_t_act_GSK3b * a['p_ApoE'] + 1) * (m_t_act_GSK3b * a['p_Ab'] + n_t_act_GSK3b) 
r_t_inactv_GSK3b = lambda a : tau_multiplier * k_t_inactv_GSK3b * a['p_GSK3b_act']
r_t_GSK3b_exp_deg = lambda a : tau_multiplier * (k_t_p_GSK3b_exp - k_t_p_GSK3b_deg * a['p_GSK3b_inact']) 
vmax_scaling_t_phos_tau = lambda a : tau_multiplier * a['p_GSK3b_act']*(a['p_GSK3b_act']/it_p_GSK3b_act)**2.547 
vmax_scaling_t_dephos_tauP = lambda a : tau_multiplier 


# ER Retraction & Collapse
r_t_LE_retro = lambda a : ER_multiplier * beta_t_LE_retro * vmax_t_LE_retro * a['p_ATP'] / (Km_t_LE_retro + a['p_ATP']) * max((mchol_t_LE_retro * a['p_chol_LE'] + nchol_t_LE_retro), 1) / dist_t_LE_trans * a['p_RTN3_axon'] * a['p_tau']/it_p_tau
r_t_LE_antero = lambda a : ER_multiplier * (vmax_t_LE_antero * a['p_ATP']) / (Km_t_LE_antero + a['p_ATP']) / dist_t_LE_trans * a['p_RTN3_PN'] * a['p_tau']/it_p_tau
r_t_RTN3_exp = lambda a : ER_multiplier * k_t_RTN3_exp
r_t_RTN3_aggregation = lambda a : ER_multiplier * 0.1 * ((a["p_RTN3_axon"] + a["p_RTN3_PN"]) - (it_p_RTN3_tot * (1 - (1-dec_t_RTN3_aggregation) * (a['p_Ab'] - it_p_Ab)/(Ab_t_RTN3_aggregation - it_p_Ab)))) 
r_t_RTN3_auto = lambda a : ER_multiplier * a["p_RTN3_HMW_cyto"] * k_t_RTN3_auto * min(1, (a["p_RTN3_axon"] / it_p_RTN3_axon)) 
r_t_RTN3_lyso = lambda a : ER_multiplier * a["p_RTN3_HMW_auto"] * k_t_RTN3_lyso * min(1, (a["p_tau"]/it_p_tau)) 
r_t_RTN3_dys_auto = lambda a : ER_multiplier * a["p_RTN3_HMW_cyto"] * k_t_RTN3_auto * max(0, (1 - a["p_RTN3_axon"] / it_p_RTN3_axon))
r_t_RTN3_dys_lyso = lambda a : ER_multiplier * a["p_RTN3_HMW_auto"] * k_t_RTN3_lyso * max(0, (1 - a["p_tau"]/it_p_tau)) 


# Calcium homeostasis rates
r_t_Ca_imp = lambda a : 1.44*1e8
r_t_mCU = lambda a : k_t_mCU1*a['p_Ca_cyto']
r_t_mNCLX = lambda a : k_t_mNCLX*a['p_Ca_mito']*(a['p_LRRK2_mut']==0) + k_t_mNCLX*0.5*a['p_Ca_mito']*(a['p_LRRK2_mut']>0)
r_t_MAM = lambda a : k_t_MAM*a['p_Ca_ER']#**1.5  #Ca/s
r_t_RyR_IP3R = lambda a : k_t_RyR_IP3R*a['p_Ca_ER']*a['p_PERK_ER']/it_p_PERK_ER 
r_t_SERCA = lambda a : k_t_SERCA_no_ATP*a['p_Ca_cyto']/a['p_NPC1_LE'] 
r_t_NCX_PMCA = lambda a: (k_t_NCX_PMCA*a['p_Ca_cyto'])*(a['p_NPC1_LE']==1) + (k_t_NCX_PMCA/5*a['p_Ca_cyto'])*(a['p_NPC1_LE']>1)


# Mitochondria

# Apoptosis / SurvivalS