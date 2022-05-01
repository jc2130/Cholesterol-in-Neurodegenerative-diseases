'''Rate functions'''

from parameters import *
from initial_tokens import *

# Cholesterol homeostasis
r_t_LDLR_endocyto = lambda a : chol_multiplier * (vmax_t_LDLR_endocyto * min([1, max([0.3,(m_t_LDLR_endocyto * (a["p_chol_ER"]/it_p_chol_ER) + n_t_LDLR_endocyto)])]) * a["p_ApoEchol_extra"]/(Km_t_LDLR_endocyto + a["p_ApoEchol_extra"]))

r_t_chol_trans_LE_ER = lambda a : chol_multiplier * k_t_chol_trans_LE_ER * a["p_chol_LE"]
r_t_chol_trans_LE_mito = lambda a : chol_multiplier * k_t_chol_trans_LE_mito * a["p_chol_LE"]
r_t_chol_trans_LE_PM = lambda a : chol_multiplier * k_t_chol_trans_LE_PM * a["p_chol_LE"]
r_t_chol_trans_ER_PM = lambda a : chol_multiplier * k_t_chol_trans_ER_PM * a["p_chol_ER"]
r_t_chol_trans_PM_ER = lambda a : chol_multiplier * k_t_chol_trans_PM_ER * a["p_chol_PM"]
r_t_chol_trans_ER_mito = lambda a : chol_multiplier * k_t_chol_trans_ER_mito * a["p_chol_ER"]  
r_t_chol_trans_PM_ECM = lambda a : chol_multiplier * k_t_chol_trans_PM_ECM * a["p_chol_PM"] * min( [3.5, max([0,(m_t_chol_trans_PM_ECM * (a["p_24OHchol_intra"]/it_p_24OHchol_intra) + n_t_chol_trans_PM_ECM)])])

r_t_27OHchol_endocyto = lambda a : chol_multiplier * disease_multiplier_27OH * k_t_27OHchol_endocyto 
r_t_24OHchol_exocyto = lambda a : chol_multiplier * k_t_24OHchol_exocyto * a["p_24OHchol_intra"] 

vmax_scaling_t_CYP27A1_metab = lambda a : chol_multiplier
vmax_scaling_t_CYP11A1_metab = lambda a : chol_multiplier
vmax_scaling_t_CYP7B1_metab = lambda a : chol_multiplier
vmax_scaling_t_CYP46A1_metab = lambda a : chol_multiplier


# ER Retraction & Collapse
r_t_LE_retro = lambda a : ER_multiplier * beta_t_LE_retro * vmax_t_LE_retro * a['p_ATP'] / (Km_t_LE_retro + a['p_ATP']) * max((mchol_t_LE_retro * a['p_chol_LE'] + nchol_t_LE_retro), 1) / dist_t_LE_trans * a['p_RTN3_axon'] * a['p_tau']/it_p_tau
r_t_LE_antero = lambda a : ER_multiplier * (vmax_t_LE_antero * a['p_ATP']) / (Km_t_LE_antero + a['p_ATP']) / dist_t_LE_trans * a['p_RTN3_PN'] * a['p_tau']/it_p_tau
r_t_RTN3_exp = lambda a : ER_multiplier * k_t_RTN3_exp

r_t_RTN3_aggregation = lambda a : ER_multiplier * 0.1 * ((a["p_RTN3_axon"] + a["p_RTN3_PN"]) - (it_p_RTN3_tot * (1 - (1-dec_t_RTN3_aggregation) * max(0,(a['p_Ab'] - it_p_Ab))/(Ab_t_RTN3_aggregation - it_p_Ab)))) # kinetics not measured and dont really matter; only equilibrium/threshold matters. just assume that 10% of excess RTN3 above threshold aggregates every second (for smoothness)

r_t_RTN3_auto = lambda a : ER_multiplier * a["p_RTN3_HMW_cyto"] * k_t_RTN3_auto * min(1, (a["p_RTN3_axon"] / it_p_RTN3_axon)) # autophagosomes per second (each containing 1 oligomer), scaled by the proportion of functional ER left
r_t_RTN3_lyso = lambda a : ER_multiplier * a["p_RTN3_HMW_auto"] * k_t_RTN3_lyso * min(1, (a["p_tau"]/it_p_tau)) # rate of transport across axon, scaled by the proportion of functional tau left
r_t_RTN3_dys_auto = lambda a : ER_multiplier * a["p_RTN3_HMW_cyto"] * k_t_RTN3_auto * max(0, (1 - a["p_RTN3_axon"] / it_p_RTN3_axon)) # any reduction in normal autophagosome formation due to lack of ER instead goes to this transition at the same rate
r_t_RTN3_dys_lyso = lambda a : ER_multiplier * a["p_RTN3_HMW_auto"] * k_t_RTN3_lyso * max(0, (1 - a["p_tau"]/it_p_tau)) # any reduction in normal lysosomal processing instead goes to this transitoin. didn't incorporate transport-inhibiting effect of dystrophic neurites yet


# Abeta Pathology 
r_t_asec_exp = lambda a : Abeta_multiplier * k_t_asec_exp * (mchol_t_asec_exp * a['p_24OHchol_intra'] + nchol_t_asec_exp) 
r_t_asec_degr = lambda a : Abeta_multiplier * k_t_asec_degr * a['p_asec'] 

r_t_APP_exp = lambda a : Abeta_multiplier * k_t_APP_exp * (1 + dis_t_APP_exp * a['p_ApoE']) #* min(max((m_t_APP_exp * a['p_ROS_mito'] + n_t_APP_exp),1),1.5)
r_t_APP_endocyto = lambda a : Abeta_multiplier * a['p_APP_pm'] * (k_t_APP_endocyto + a['p_ApoE'] * dis_t_APP_endocyto)
r_t_APP_endo_event = lambda a :Abeta_multiplier *  k_t_APP_endo_event * a['p_APP_endo']

r_t_bsec_exp = lambda a : Abeta_multiplier * k_t_bsec_exp * (mchol_t_bsec_exp * a['p_27OHchol_intra'] + nchol_t_bsec_exp) * (mRTN_t_bsec_exp * a['p_RTN3_axon'] + nRTN_t_bsec_exp) #* max(min((mROS_t_bsec_exp * a['p_ROS_mito'] + nROS_t_bsec_exp),2),1) #arbitrary upper limit of 2 to avoid breaking model in extreme runs
r_t_bsec_degr = lambda a : Abeta_multiplier * k_t_bsec_degr * a['p_bsec'] 

r_t_gsec_exp = lambda a : Abeta_multiplier * k_t_gsec_exp # need to incorporate effect of p_ROS_mito
r_t_gsec_degr = lambda a : Abeta_multiplier * k_t_gsec_degr * a['p_gsec'] 

r_t_Ab_degr = lambda a : Abeta_multiplier * k_t_Ab_degr * a['p_Ab']

vmax_scaling_t_APP_asec_cleav = lambda a : Abeta_multiplier * a['p_asec']
vmax_scaling_t_APP_bsec_cleav = lambda a : Abeta_multiplier * a['p_bsec'] * (mchol_t_APP_bsec_cleav * a['p_chol_PM'] + nchol_t_APP_bsec_cleav) * (1 + a['p_age'] * age_t_APP_bsec_cleav)
vmax_scaling_t_CTF99_gsec_cleav = lambda a : Abeta_multiplier * a['p_gsec']


# Tau Pathology
r_t_actv_GSK3b = lambda a : tau_multiplier * k_t_actv_GSK3b * a['p_GSK3b_inact'] * (dis_t_act_GSK3b * a['p_ApoE'] + 1) * (m_t_act_GSK3b * a['p_Ab'] + n_t_act_GSK3b) 
r_t_inactv_GSK3b = lambda a : tau_multiplier * k_t_inactv_GSK3b * a['p_GSK3b_act']
r_t_GSK3b_exp_deg = lambda a : tau_multiplier * (k_t_p_GSK3b_exp - k_t_p_GSK3b_deg * a['p_GSK3b_inact']) 

vmax_scaling_t_phos_tau = lambda a : tau_multiplier * a['p_GSK3b_act']*(a['p_GSK3b_act']/it_p_GSK3b_act)**2.547 #(1*a['p_cas']+1) * (0.2*a['p_GSK3b_act']+1))
vmax_scaling_t_dephos_tauP = lambda a : tau_multiplier 


# Calcium Homeostasis
r_t_Ca_imp = lambda a : 1.44*1e8
r_t_mCU = lambda a : k_t_mCU1 *a['p_Ca_cyto'] #*(a['p_Ca_mito']<=3*1e7) + (a['p_Ca_mito']>3*1e7)*k_t_mCU2*a['p_Ca_cyto']
r_t_mNCLX = lambda a : k_t_mNCLX*a['p_Ca_mito'] #Ca/s #margaux: I removed the dependancy on lRRK2 compared to what PD had done
r_t_MAM = lambda a : k_t_MAM*a['p_Ca_ER'] #Ca/s #/(1+(a['p_Ca_mito']-it_p_Ca_mito))
r_t_RyR_IP3R = lambda a : k_t_RyR_IP3R*a['p_Ca_ER'] #100 #Ca/s
#WITHOUT ATP
r_t_SERCA = lambda a : k_t_SERCA_no_ATP*a['p_Ca_cyto']#Ca/s
#WITH ATP
#r_t_SERCA = lambda a : k_t_SERCA*a['p_Ca_cyto']*a['p_ATP'] #Ca/s
r_t_NCX_PMCA = lambda a: (k_t_NCX_PMCA*a['p_Ca_cyto'])*(a['p_Ca_cyto']<=7.5*1e7)+ (k_t_NCX_PMCA/10*a['p_Ca_cyto'])*(a['p_Ca_cyto']>7.5*1e7)


# Energy metabolism
r_t_krebs = lambda a : k_t_krebs * a['p_ADP'] * a['p_Ca_mito'] * min( [1.0, max([0.0018,(m_t_ETC_inhib_Ab * a["p_Ab"] + n_t_ETC_inhib_Ab)])]) 
r_t_ATP_hydro_mito = lambda a : k_t_ATP_hydro_mito * a['p_ATP']
r_t_ETC = lambda a : k_t_ETC * a['p_ADP'] * a['p_Ca_mito'] *  a['p_reduc_mito'] / a['p_ROS_mito']**0.5 / max(a['p_chol_mito'], it_p_chol_mito) * min( [1.0, max([0.0018,(m_t_ETC_inhib_Ab * a["p_Ab"] + n_t_ETC_inhib_Ab)])])
r_t_ROS_metab = lambda a : k_t_ROS_metab * a['p_ROS_mito'] / max(a['p_chol_mito'], it_p_chol_mito)
r_t_mito_dysfunc = lambda a : k_t_mito_dysfunc * (m_t_mito_dysfunc * a["p_ROS_mito"] + n_t_mito_dysfunc) * (m_t_mito_dysfunc_Ab * a["p_Ab"] + n_t_mito_dysfunc_Ab)
r_t_cas3_inact = lambda a : k_t_cas3_inact * a["p_cas3"]
r_t_ROS_gener_Ab = lambda a : k_t_ROS_gener_Ab * (a["p_Ab"] - it_p_Ab)