'''Initial tokens'''

# Cholesterol homeostasis
it_p_ApoEchol_extra = 1.60e6
it_p_ApoEchol_EE = 2.5 * 10 **3 
it_p_chol_LE = 2.65 * 10 ** 9 
it_p_chol_mito = 3.79*10**9
it_p_chol_ER = 4.37*10**8
it_p_chol_PM = 8.45*10**10
it_p_24OHchol_extra = 4.65*10**3
it_p_24OHchol_intra = 3.98*10**8
it_p_27OHchol_intra = 3.06*10**7
it_p_27OHchol_extra = 1.72*10**3
it_p_7HOCA = 0
it_p_preg = 0
it_p_SNCA_act_extra = 113    # PD specific place
it_p_SNCAApoEchol_extra = 0  # PD specific place
it_p_SNCAApoEchol_intra = 0  # PD specific place
it_p_NPC1_LE = 1 # 1 for healthy, 660000 for mutated
it_p_PERK_ER = 1*10**6 * ((it_p_NPC1_LE < 5) * 1 + (it_p_NPC1_LE > 5) * 100)
it_p_TRADD_cyto = 298540


# Lewy body pathology
it_p_SNCA_act = 4.3383e8
it_p_SNCA_inact = 0
it_p_SNCA_olig = 0
it_p_LB = 0
it_p_Fe2 = 2e8


# PD specific Mutations (discrete places)
it_p_GBA1 = 0
it_p_ABCGA1 = 0
it_p_LRRK2_mut = 0
it_p_VPS35 = 0
it_p_DJ1 = 0
it_p_LAMP2A = 0 #therapeutic PD, set to 1 for PD drug


# ER Retraction & Collapse
it_p_RTN3_axon =   5344426
it_p_RTN3_PN =   6655574
it_p_RTN3_tot = it_p_RTN3_axon + it_p_RTN3_PN # not actually for initializing a place, just useful in some calculations
it_p_RTN3_HMW_cyto = 0
it_p_RTN3_HMW_auto = 0
it_p_RTN3_HMW_lyso = 0
it_p_RTN3_HMW_dys1 = 0
it_p_RTN3_HMW_dys2 = 0


# Abeta Pathology
it_p_asec = 6.05e6
it_p_APP_pm = 3.49e6
it_p_sAPPa = 0 # sink
it_p_CTF83 = 0 # sink
it_p_APP_endo = 1.91e6
it_p_bsec = 6.73e5
it_p_sAPPb = 0 # sink
it_p_CTF99 = 1930 #1.68E6 
it_p_gsec = 3.37e6
it_p_AICD = 0 # sink
it_p_Ab = 34048 # 
it_p_ApoE = 0 # 0 if healthy/non-risk genotype, 1 for ApoE4 risk factor allele
it_p_age = 0 # 0 if young/no age risk factor, 1 representing 80 years of age risk factor (affects BACE1 activity)


# Tau Pathology
it_p_GSK3b_inact = 2684994
it_p_GSK3b_act = 2815005
it_p_tauP = 6.36e6
it_p_tau = 3.72e7
it_p_PP_2A = 1.33e7


# Calcium Homeostasis
it_p_Ca_cyto = 7.2*1e5
it_p_Ca_mito = 3e7 
it_p_Ca_ER = 1.8*1e9
it_p_ADP = 5.42*10**6
it_p_ATP = 5.42*10**9 #1*1e14 
it_p_Ca_LE = 2.44*10**5

# Energy Metabolism
it_p_cas3 = 13014.7584
it_p_ATP = 5.42*10**9
it_p_ADP = 5.42*10**6
it_p_reduc_mito = 2.71*10**9
it_p_ROS_mito = 13557.04
it_p_H2O_mito = 0

#NPC Therapeutic
it_p_NAC = 0