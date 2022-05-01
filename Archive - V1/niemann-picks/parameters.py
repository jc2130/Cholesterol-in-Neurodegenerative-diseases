# Rate constants for all transitions
from initial_tokens import *

# multiplicative rate factors for increasing rates of slow modules
ER_multiplier = 0.5
Chol_speed_up_parameter = 300


# Choleserol Homeostasis
m_t_LDLR_endocyto = -1.0682
n_t_LDLR_endocyto = 2.0682
Km_t_LDLR_endocyto = 1.30 * 10 ** 6
vmax_t_LDLR_endocyto = Chol_speed_up_parameter * 3.61633 * 10 ** 4
st_t_CYP27A1_metab = 0.13158 # CYP27A1-accessible portion of mitochondrial cholesterol (to scale Km)
Km_t_CYP27A1_metab = 4.77 * 10 ** 7 / st_t_CYP27A1_metab
vmax_t_CYP27A1_metab = Chol_speed_up_parameter * 2.56 * 10 ** 3
st_t_CYP11A1_metab = 0.13158 # CYP11A1-accessible portion of mitochondrial cholesterol (to scale Km)
Km_t_CYP11A1_metab = 7.59 * 10 ** 7 / st_t_CYP11A1_metab 
vmax_t_CYP11A1_metab = Chol_speed_up_parameter * 6.35 * 10 ** 4
k_t_chol_trans_LE_ER = Chol_speed_up_parameter * 2.55357 * 10 ** (-4) 
k_t_chol_trans_LE_mito = Chol_speed_up_parameter * 2.36 * 10 ** (-9) 
k_t_chol_trans_LE_PM = Chol_speed_up_parameter * 0.002406761 
k_t_chol_trans_ER_PM = Chol_speed_up_parameter * 1.725 * 10 ** (-3) 
k_t_chol_trans_PM_ER = Chol_speed_up_parameter * 1.56 * 10 ** (-6) 
k_t_chol_trans_ER_mito = Chol_speed_up_parameter * 1.1713 * 10 ** (-7) 
k_t_chol_exocyto = Chol_speed_up_parameter * 8.2859 * 10 ** (-5) 

# Tau Pathology
k_t_phos_tau = 5.07488*10**(-18) 
k_t_dephos_tauP = 0.004033267867924528301886792452830188679245283018867924528
k_t_actv_GSK3b = 8.33e-3
m_t_act_GSK3b = 4.07e-7 
n_t_act_GSK3b = 1 - m_t_act_GSK3b * it_p_Ab
dis_t_act_GSK3b = 0.433
k_t_inactv_GSK3b = 7.95e-3
Km_t_phos_tau = 9.22e7
kcat_t_phos_tau = 0.146464095 
Km_t_dephos_tauP = 6.29e7
vmax_t_dephos_tauP = 1.17*1.1e6  
k_t_p_GSK3b_deg = 100*1.6e-5 
k_t_p_GSK3b_exp = k_t_p_GSK3b_deg * it_p_GSK3b_inact


# ER Retraction & Collapse
beta_t_LE_retro = 1.667 
dist_t_LE_trans = 75e4 
mchol_t_LE_retro = 2.27e-9 
nchol_t_LE_retro = 1 - mchol_t_LE_retro * it_p_chol_LE 
vmax_t_LE_retro = 892 
Km_t_LE_retro = 3510864000 
vmax_t_LE_antero = 814 
Km_t_LE_antero = 614040000 
ATPcons_t_LE_trans = 0 
k_t_RTN3_exp = 113.3
Ab_t_RTN3_aggregation = 641020
dec_t_RTN3_aggregation = 0.762
k_t_RTN3_auto = 0.011111111
k_t_RTN3_lyso = 0.000826667
mitprop_t_RTN3_dys_auto = 0.885


# Calcium and Sphingosine Homeostasis
k_t_RyR_IP3R = 100/(1.8*1e9)
k_t_Ca_cyto_LE = 0.5*10**(-2)
k_t_TRPML1 = 0.92*10**1
k_t_Ca_imp = 1.44*10**7
k_t_NCX_PMCA = 10
k_t_mNCLX = 0.066666667
k_t_SERCA_no_ATP=0.05638 
k_t_SERCA=k_t_SERCA_no_ATP/5.42e9 
k_t_mCU1=(1*1e6)/(17854326)
k_t_MAM = (1*1e6)/(1.8*1e9)
k_t_NaK_ATPase = 7.28448*10**8*2


# Mitochondria
k_t_krebs = 1.32*10**(-8)                                    
k_t_ATP_hydro_mito = 1.92*10**(-1)*0.3
k_t_ETC = 4.419509028*10**(-5)                
k_t_ROS_metab = 3.078148*10**(15)                        


# Apoptosis / Survival
k_t_ERstress_Ca_cyto_induc = 1.582602*10**(-6)
k_t_BAX_actv = 1
k_t_Bcl2_actv = 360000
k_t_mito_dysfunc = 3380.67
k_t_mito_dysfunc_2 = 1
k_t_TRADD_actv = 355.9059759501279
k_t_Cas3_actv = 304.5 / 27140
k_t_Cas3_actv_2 = 303.35 / 27140
k_t_Cas3_degr = 27325.972
k_t_mTORC1_actv = 13.2
k_t_NFkB_actv = 1.713*10**(-7)
k_t_NFkB_actv_2 = 4.051*10**(-6)
k_t_NFkB_actv_3 = 1.485*10**(-5)
k_t_NFkB_degr = 2.7572*10**2                          