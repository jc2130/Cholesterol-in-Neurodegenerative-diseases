# NPCD HFPN subnets


# Imports
import os
import sys
import numpy as np

cwd = os.getcwd()

root_folder = os.sep+"team-project"
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep+"utils"+os.sep)

# Import class to work with hybrid functional Petri nets
from hfpn import HFPN

# Import parameters, rate equations, initial tokens and firing conditions
from parameters import *
from rate_functions import *
from initial_tokens import *
from firing_conditions import *


### Cholesterol Homeostasis
def add_cholesterol_homeostasis(hfpn):

    ## Places
    hfpn.add_place(initial_tokens=it_p_ApoEchol_extra, place_id="p_ApoEchol_extra", label="Conc of chol-ApoE being endocytosized", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_mito, place_id="p_chol_mito", label="Chol conc in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_NPC1_LE, place_id="p_NPC1_LE", label="Mutated NPC1 in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_ER, place_id="p_chol_ER", label="Chol conc in ER", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_PM, place_id="p_chol_PM", label="cholesterol in PM", continuous=True)
    hfpn.add_place(initial_tokens=it_p_preg, place_id="p_preg", label="Pregnenolon", continuous=True)
    hfpn.add_place(initial_tokens=it_p_SINK, place_id="p_SINK", label="Sink for cholesterol export from PM", continuous=True)
    #There is no real needed output place, so we just create artificial "Sink" for flow of cholestrol out of PM

    #Places from other modules
    hfpn.add_place(initial_tokens=it_p_PERK_ER, place_id="p_PERK_ER", label="Conc of PERK in ER", continuous=True)


    ## Transitions
    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_LE_PM",
                label                          = "Chol transport LE-PM",
                input_place_ids                 = ["p_chol_LE"],
                firing_condition             = fc_t_chol_trans_LE_PM, 
                reaction_speed_function         = r_t_chol_trans_LE_PM,
                consumption_coefficients     = [1],
                output_place_ids             = ["p_chol_PM"],
                production_coefficients         = [1])

    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_ER_PM",
                label                          = "Chol transport ER-PM",
                input_place_ids                 = ["p_chol_ER"],
                firing_condition             = fc_t_chol_trans_ER_PM,
                reaction_speed_function         = r_t_chol_trans_ER_PM,
                consumption_coefficients     = [1],
                output_place_ids             = ["p_chol_PM"],
                production_coefficients         = [1])

    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_PM_ER",
                label                          = "Chol transport PM-ER",
                input_place_ids                 = ["p_chol_PM"],
                firing_condition             = fc_t_chol_trans_PM_ER,
                reaction_speed_function         = r_t_chol_trans_PM_ER,
                consumption_coefficients     = [1],
                output_place_ids             = ["p_chol_ER"],
                production_coefficients         = [1])

    hfpn.add_transition(  transition_id = 't_chol_exocyto', 
                    label = 'Exocytosis of cholesterol', 
                    input_place_ids = ['p_chol_PM'],
                    firing_condition = fc_t_chol_exocyto,
                    consumption_speed_functions = [lambda a : k_t_chol_exocyto * a["p_chol_PM"]],
                    output_place_ids = ['p_SINK'],
                    production_speed_functions = [lambda a : 0])
                    

    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_LDLR_endocyto",
                label                          = "LDLR endocyto",
                input_place_ids                 = ["p_ApoEchol_extra", "p_chol_ER"],
                firing_condition             = fc_t_LDLR_endocyto,
                reaction_speed_function         = r_t_LDLR_endocyto, 
                consumption_coefficients     = [0,0],
                output_place_ids             = ["p_chol_LE"],
                production_coefficients         = [354])

    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_LE_ER",
                label                          = "Chol transport LE-ER",
                input_place_ids                 = ["p_chol_LE", "p_NPC1_LE"],
                firing_condition             = fc_t_chol_trans_LE_ER,
                reaction_speed_function         = r_t_chol_trans_LE_ER,
                consumption_coefficients     = [1,0],
                output_place_ids             = ["p_chol_ER"],
                production_coefficients         = [1])

    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_LE_mito",
                label                          = "Chol transport LE-mito",
                input_place_ids                 = ["p_chol_LE", "p_PERK_ER"],
                firing_condition             = fc_t_chol_trans_LE_mito,
                reaction_speed_function         = r_t_chol_trans_LE_mito,
                consumption_coefficients     = [1,0], 
                output_place_ids             = ["p_chol_mito"],
                production_coefficients         = [1])


    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_chol_trans_ER_mito",
                label                          = "Chol transport ER-mito",
                input_place_ids                 = ["p_chol_ER", "p_PERK_ER"],
                firing_condition             = fc_t_chol_trans_ER_mito,
                reaction_speed_function         = r_t_chol_trans_ER_mito,
                consumption_coefficients     = [1,0],
                output_place_ids             = ["p_chol_mito"],
                production_coefficients         = [1])


    hfpn.add_transition_with_speed_function(
                transition_id                 = "t_CYP11A1_metab",
                label                          = "Chol metab CYP11A1",
                input_place_ids                 = ["p_chol_mito"],
                firing_condition             = fc_t_CYP11A1_metab,
                reaction_speed_function         = r_t_CYP11A1_metab,
                consumption_coefficients     = [1],
                output_place_ids             = ["p_preg"],
                production_coefficients         = [1])

### Tau Pathology
def add_tau_pathology(hfpn):

    ## Places
    hfpn.add_place(initial_tokens=it_p_tau, place_id="p_tau", label="Tau in cyto", continuous=True)

    hfpn.add_place(initial_tokens=it_p_tauP, place_id="p_tauP", label="Tau phosphorylated in cyto", continuous=True)

    hfpn.add_place(initial_tokens=it_p_GSK3b_inact, place_id='p_GSK3b_inact', label='Inactive GSK3 Beta kinase', continuous=True)

    hfpn.add_place(initial_tokens=it_p_GSK3b_act, place_id='p_GSK3b_act', label='Active GSK3 Beta kinase', continuous=True)

    #Places from other modules
    hfpn.add_place(initial_tokens=it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol in LE", continuous=True)

    ## Transitions
    hfpn.add_transition(  transition_id = 't_phos_tau', 
                    label = 'phosphorylation of tau', 
                    input_place_ids = ['p_chol_LE', 'p_tau'],
                    firing_condition = fc_t_phos_tau,
                    consumption_speed_functions = [lambda a : 0,
                                                   lambda a: k_t_phos_tau * a['p_tau'] * a['p_chol_LE']**(1.5)],
                    output_place_ids = ['p_tauP'],  
                    production_speed_functions = [lambda a : k_t_phos_tau * a['p_tau'] * a['p_chol_LE']**(1.5)])

    hfpn.add_transition(  transition_id = 't_dephos_tauP', 
                    label = 'dephosphorylation of tau', 
                    input_place_ids = ['p_tauP'],
                    firing_condition = fc_t_dephos_tauP,
                    consumption_speed_functions = [lambda a: k_t_dephos_tauP * a['p_tauP']],
                    output_place_ids = ['p_tau'],  
                    production_speed_functions = [lambda a : k_t_dephos_tauP * a['p_tauP']])


### ER Retraction & Collapse
def add_ER_retraction_collapse(hfpn):

    ## Places 
    hfpn.add_place(it_p_RTN3_axon, place_id="p_RTN3_axon", label="monomeric RTN3 (axonal)", continuous=True)

    hfpn.add_place(it_p_RTN3_PN, place_id="p_RTN3_PN", label="monomeric RTN3 (perinuclear)", continuous=True)

    hfpn.add_place(it_p_RTN3_HMW_cyto, place_id="p_RTN3_HMW_cyto", label="HMW RTN3 (cytosol)", continuous=True)

    hfpn.add_place(it_p_RTN3_HMW_auto, place_id="p_RTN3_HMW_auto", label="HMW RTN3 (autophagosome)", continuous=True)

    hfpn.add_place(it_p_RTN3_HMW_lyso, place_id="p_RTN3_HMW_lyso", label="HMW RTN3 (degraded in lysosome)", continuous=True)

    hfpn.add_place(it_p_RTN3_HMW_dys1, place_id="p_RTN3_HMW_dys1", label="HMW RTN3 (type I/III dystrophic neurites)", continuous=True)

    hfpn.add_place(it_p_RTN3_HMW_dys2, place_id="p_RTN3_HMW_dys2", label="HMW RTN3 (type II dystrophic neurites)", continuous=True)

    # From other modules
    hfpn.add_place(initial_tokens=it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_tau, place_id="p_tau", label="Tau in cyto", continuous=True)
    hfpn.add_place(initial_tokens=it_p_ATP, place_id="p_ATP", label="Conc of ATP in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_ADP, place_id="p_ADP", label="Conc of ADP in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_Ab, place_id='p_Ab', label='', continuous=True)


    ## Transitions
    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_exp',
                        label = 'Expression rate of RTN3',
                        input_place_ids = [], 
                        firing_condition = fc_t_RTN3_exp,
                        reaction_speed_function = r_t_RTN3_exp, 
                        consumption_coefficients = [],
                        output_place_ids = ['p_RTN3_PN'],
                        production_coefficients = [1])
    
    hfpn.add_transition_with_speed_function(
                        transition_id = 't_LE_retro',
                        label = 'retrograde transport of LEs & ER',
                        input_place_ids = ['p_ATP','p_chol_LE','p_RTN3_axon', 'p_tau'], 
                        firing_condition = fc_t_LE_retro,
                        reaction_speed_function = r_t_LE_retro, 
                        consumption_coefficients = [ATPcons_t_LE_trans, 0, 1, 0], 
                        output_place_ids = ['p_ADP','p_RTN3_PN'],
                        production_coefficients = [ATPcons_t_LE_trans, 1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_LE_antero',
                        label = 'anterograde transport of LEs & ER',
                        input_place_ids = ['p_ATP','p_RTN3_PN', 'p_tau'], 
                        firing_condition = fc_t_LE_antero,
                        reaction_speed_function = r_t_LE_antero, 
                        consumption_coefficients = [ATPcons_t_LE_trans, 1, 0], 
                        output_place_ids = ['p_ADP','p_RTN3_axon'],
                        production_coefficients = [ATPcons_t_LE_trans, 1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_aggregation',
                        label = 'aggregation of monomeric RTN3 into HMW RTN3',
                        input_place_ids = ['p_RTN3_axon', 'p_RTN3_PN', 'p_Ab'], 
                        firing_condition = fc_t_RTN3_aggregation, 
                        reaction_speed_function = r_t_RTN3_aggregation,
                        consumption_coefficients = [1, 1, 0],
                        output_place_ids = ['p_RTN3_HMW_cyto'],
                        production_coefficients = [1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_auto',
                        label = 'functional autophagy of HMW RTN3',
                        input_place_ids = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'], 
                        firing_condition = fc_t_RTN3_auto, 
                        reaction_speed_function = r_t_RTN3_auto,
                        consumption_coefficients = [1, 0],
                        output_place_ids = ['p_RTN3_HMW_auto'],
                        production_coefficients = [1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_lyso',
                        label = 'functional delivery of HMW RTN3 to the lysosome',
                        input_place_ids = ['p_RTN3_HMW_auto', 'p_tau'], 
                        firing_condition = fc_t_RTN3_lyso, 
                        reaction_speed_function = r_t_RTN3_lyso,
                        consumption_coefficients = [1, 0],
                        output_place_ids = ['p_RTN3_HMW_lyso'],
                        production_coefficients = [1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_dys_auto',
                        label = 'dysfunctional autophagy of HMW RTN3',
                        input_place_ids = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'], 
                        firing_condition = fc_t_RTN3_dys_auto, 
                        reaction_speed_function = r_t_RTN3_dys_auto,
                        consumption_coefficients = [1, 0],
                        output_place_ids = ['p_RTN3_HMW_dys1'],
                        production_coefficients = [1]) 

    hfpn.add_transition_with_speed_function(
                        transition_id = 't_RTN3_dys_lyso',
                        label = 'dysfunctional delivery of HMW RTN3 to the lysosome',
                        input_place_ids = ['p_RTN3_HMW_auto', 'p_RTN3_HMW_dys1', 'p_tau'], 
                        firing_condition = fc_t_RTN3_dys_lyso, 
                        reaction_speed_function = r_t_RTN3_dys_lyso,
                        consumption_coefficients = [1, 0, 0],
                        output_place_ids = ['p_RTN3_HMW_dys2'],
                        production_coefficients = [1]) 
##Calcium
def add_calcium_homeostasis(hfpn):

    # Add places
    hfpn.add_place(initial_tokens=it_p_Ca_LE, place_id="p_Ca_LE", label="Ca conc in LE", continuous=True)

    hfpn.add_place(initial_tokens=it_p_Ca_cyto, place_id="p_Ca_cyto", label="Ca conc in cyto", continuous=True)

    hfpn.add_place(initial_tokens=it_p_Ca_ER, place_id="p_Ca_ER", label="Ca conc in ER", continuous=True)

    hfpn.add_place(initial_tokens=it_p_Ca_mito, place_id="p_Ca_mito", label="Ca conc in mito", continuous=True)

    #ON/OFF switches for pacemaking:
    hfpn.add_place(initial_tokens=it_p_Ca_extra, place_id="p_Ca_extra", label="Extracell. Ca conc imported via L-channel", continuous=False)

    hfpn.add_place(initial_tokens=0, place_id="p_on2", label="on2", continuous=False)

    hfpn.add_place(initial_tokens=0, place_id="p_on3", label="on3", continuous=False)

    hfpn.add_place(initial_tokens=0, place_id="p_on4", label="on4", continuous=False)

    #PD specific, remove after integration
    hfpn.add_place(it_p_LRRK2_mut, "p_LRRK2_mut","LRRK2 - mutated", continuous = True)

    # Places of other modules
    hfpn.add_place(initial_tokens=it_p_PERK_ER, place_id="p_PERK_ER", label="Conc of PERK in ER", continuous=True)
    hfpn.add_place(initial_tokens=it_p_NPC1_LE, place_id="p_NPC1_LE", label="Mutated NPC1 in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_ATP, place_id="p_ATP", label="Conc of ATP in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_ADP, place_id="p_ADP", label="Conc of ADP in mito", continuous=True)


    ## Transition
    hfpn.add_transition_with_speed_function(  
                    transition_id = 't_RyR_IP3R',
                    label = 'Ca export from ER',
                    input_place_ids = ['p_Ca_extra','p_Ca_ER','p_PERK_ER'],
                    firing_condition = fc_t_RyR_IP3R,
                    reaction_speed_function = r_t_RyR_IP3R,
                    consumption_coefficients = [0,1,0], 
                    output_place_ids = ['p_Ca_cyto'],         
                    production_coefficients = [1]) 

    hfpn.add_transition(  transition_id = 't_Ca_cyto_LE', 
                    label = 'NAADP facilitated transport of Ca from cyto to LE', 
                    input_place_ids = ['p_NPC1_LE', 'p_Ca_cyto'],
                    firing_condition = fc_t_Ca_cyto_LE,
                    consumption_speed_functions = [lambda a : 0,
                                                   lambda a : k_t_Ca_cyto_LE * a['p_Ca_cyto'] / a['p_NPC1_LE']**(0.2)],
                    output_place_ids = ['p_Ca_LE'],  
                    production_speed_functions = [lambda a : k_t_Ca_cyto_LE * a['p_Ca_cyto'] / a['p_NPC1_LE']**(0.2)])

    hfpn.add_transition(  transition_id = 't_TRPML1', 
                    label = 'TRPML1 facilitated transport of Ca from LE to cyto', 
                    input_place_ids = ['p_Ca_LE', 'p_chol_LE'],
                    firing_condition = fc_t_TRPML1,
                    consumption_speed_functions = [lambda a : k_t_TRPML1*a['p_Ca_LE']/a['p_chol_LE']**(0.154),
                                                   lambda a : 0],
                    output_place_ids = ['p_Ca_cyto'],  
                    production_speed_functions = [lambda a : k_t_TRPML1*a['p_Ca_LE']/a['p_chol_LE']**(0.154)])

    hfpn.add_transition_with_speed_function(  
                    transition_id = 't_Ca_imp',
                    label = 'L-type Ca channel',
                    input_place_ids = ['p_Ca_extra'],
                    firing_condition = fc_t_Ca_imp,
                    reaction_speed_function = r_t_Ca_imp,
                    consumption_coefficients = [0], 
                    output_place_ids = ['p_Ca_cyto'],         
                    production_coefficients = [1])  

    hfpn.add_transition_with_speed_function(  
                    transition_id = 't_NCX_PMCA',
                    label = 'Ca efflux to extracellular space',
                    input_place_ids = ['p_Ca_cyto','p_on3','p_NPC1_LE'],
                    firing_condition = fc_t_NCX_PMCA,
                    reaction_speed_function = r_t_NCX_PMCA,
                    consumption_coefficients = [1,0,0], 
                    output_place_ids = [],       
                    production_coefficients = []) 

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_mNCLX',
                    label = 'Ca export from mitochondria via mNCLX',
                    input_place_ids = ['p_Ca_mito','p_LRRK2_mut'],
                    firing_condition = fc_t_mNCLX,
                    reaction_speed_function = r_t_mNCLX,
                    consumption_coefficients = [1,0], 
                    output_place_ids = ['p_Ca_cyto'],         
                    production_coefficients = [1])

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_SERCA',
                    label = 'Ca import to ER',
                    input_place_ids = ['p_Ca_cyto','p_ATP','p_NPC1_LE'],
                    firing_condition = fc_t_SERCA,
                    reaction_speed_function = r_t_SERCA,
                    consumption_coefficients = [1,1,0], 
                    output_place_ids = ['p_Ca_ER','p_ADP'],         
                    production_coefficients = [1,1]) 

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_mCU',
                    label = 'Ca import into mitochondria via mCU',
                    input_place_ids = ['p_Ca_cyto','p_Ca_mito'],
                    firing_condition = fc_t_mCU,
                    reaction_speed_function = r_t_mCU,
                    consumption_coefficients = [1,0], 
                    output_place_ids = ['p_Ca_mito'],         
                    production_coefficients = [1]) 

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_MAM',
                    label = 'Ca transport from ER to mitochondria',
                    input_place_ids = ['p_Ca_ER','p_Ca_mito'],
                    firing_condition = fc_t_MAM,
                    reaction_speed_function = r_t_MAM,
                    consumption_coefficients = [1,0], 
                    output_place_ids = ['p_Ca_mito'],         
                    production_coefficients = [1]) 

    #ON/OFF switches for pacemaking:
    hfpn.add_transition_with_speed_function(
                    transition_id = 't_A',
                    label = 'A',
                    input_place_ids = ['p_on4'],
                    firing_condition = lambda a: a['p_on4']==1,
                    reaction_speed_function = lambda a: 1,
                    consumption_coefficients = [1], 
                    output_place_ids = ['p_Ca_extra'],         
                    production_coefficients = [1],
                    delay=0.5)

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_B',
                    label = 'B',
                    input_place_ids = ['p_Ca_extra'],
                    firing_condition = lambda a: a['p_Ca_extra']==1,
                    reaction_speed_function = lambda a: 1,
                    consumption_coefficients = [1], 
                    output_place_ids = ['p_on2'],         
                    production_coefficients = [1],
                    delay=0.5) 

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_C',
                    label = 'C',
                    input_place_ids = ['p_on2'],
                    firing_condition = lambda a: a['p_on2']==1,
                    reaction_speed_function = lambda a: 1,
                    consumption_coefficients = [1], 
                    output_place_ids = ['p_on3'],         
                    production_coefficients = [1],
                    delay=0) 

    hfpn.add_transition_with_speed_function(
                    transition_id = 't_D',
                    label = 'D',
                    input_place_ids = ['p_on3'],
                    firing_condition = lambda a: a['p_on3']==1,
                    reaction_speed_function = lambda a: 1,
                    consumption_coefficients = [1], 
                    output_place_ids = ['p_on4'],         
                    production_coefficients = [1],
                    delay=0.5)

    # Link to energy metabolism in that it needs ATP replenishment
    hfpn.add_transition(  transition_id = 't_NaK_ATPase', 
                    label = 'NaK ATPase', 
                    input_place_ids = ['p_ATP', 'p_on3'],
                    firing_condition = fc_t_NaK_ATPase,
                    consumption_speed_functions = [lambda a : k_t_NaK_ATPase,
                                                   lambda a : 0],
                    output_place_ids = ['p_ADP'],  
                    production_speed_functions = [lambda a : k_t_NaK_ATPase])


##Mito
def add_mitochondria(hfpn):

    ## Places
    hfpn.add_place(initial_tokens=it_p_ATP, place_id="p_ATP", label="Conc of ATP in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_ADP, place_id="p_ADP", label="Conc of ADP in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_reduc_mito, place_id="p_reduc_mito", label="Conc of reducing agents (NADH, FADH) in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_ROS_mito, place_id="p_ROS_mito", label="Conc of ROS in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_H2O_mito, place_id="p_H2O_mito", label="Conc of H2O in mito", continuous=True)

    # Places from other modules
    hfpn.add_place(initial_tokens=it_p_Ca_mito, place_id="p_Ca_mito", label="Ca conc in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_mito, place_id="p_chol_mito", label="Chol conc in mito", continuous=True)

    ## Transitions
    hfpn.add_transition(  transition_id = 't_krebs', 
                    label = 'Krebs cycle', 
                    input_place_ids = ['p_ADP', 'p_Ca_mito'],
                    firing_condition = fc_t_krebs,
                    consumption_speed_functions = [lambda a : k_t_krebs * a['p_ADP'] * a['p_Ca_mito'],
                                                   lambda a : 0],
                    output_place_ids = ['p_reduc_mito', 'p_ATP'],  
                    production_speed_functions = [lambda a : k_t_krebs * a['p_ADP'] * a['p_Ca_mito'] * 4,
                                                  lambda a : k_t_krebs * a['p_ADP'] * a['p_Ca_mito']])

    hfpn.add_transition(  transition_id = 't_ATP_hydro_mito', 
                    label = 'ATP hydrolysis by cellular processes', 
                    input_place_ids = ['p_ATP'],
                    firing_condition = fc_t_ATP_hydro_mito,
                    consumption_speed_functions = [lambda a : k_t_ATP_hydro_mito * a['p_ATP']], 
                    output_place_ids = ['p_ADP'],  
                    production_speed_functions = [lambda a : k_t_ATP_hydro_mito * a['p_ATP']])

    hfpn.add_transition(  transition_id = 't_ETC', 
                    label = 'Electron transport chain', 
                    input_place_ids = ['p_reduc_mito', 'p_ADP', 'p_Ca_mito', 'p_chol_mito', 'p_ROS_mito'], 
                    firing_condition = fc_t_ETC,
                    consumption_speed_functions = [lambda a : 6/242 * 22 / 3 * k_t_ETC * a['p_ADP'] * a['p_Ca_mito'] *  a['p_reduc_mito'] / a['p_chol_mito'] / a['p_ROS_mito']**0.5,
                                                   lambda a : 22 * k_t_ETC * a['p_ADP'] * a['p_Ca_mito'] * a['p_reduc_mito'] / a['p_chol_mito']  / a['p_ROS_mito']**0.5,
                                                   lambda a : 0,
                                                   lambda a : 0,
                                                   lambda a : 0],
                    output_place_ids = ['p_ATP', 'p_ROS_mito'],  
                    production_speed_functions = [lambda a : 22 * k_t_ETC * a['p_ADP'] * a['p_Ca_mito'] * a['p_reduc_mito'] / a['p_chol_mito']  / a['p_ROS_mito']**0.5,
                                                  lambda a : 0.02 * 2 / 2 * 6 / 242 * 22 / 3 * k_t_ETC * a['p_ADP'] * a['p_Ca_mito'] * a['p_reduc_mito'] / a['p_chol_mito']  / a['p_ROS_mito']**0.5])

    hfpn.add_transition(  transition_id = 't_ROS_metab', 
                    label = 'Oxidation of proteins/lipids involved in ETC', 
                    input_place_ids = ['p_ROS_mito', 'p_chol_mito'],
                    firing_condition = fc_t_ROS_metab,
                    consumption_speed_functions = [lambda a : k_t_ROS_metab * a['p_ROS_mito'] / a['p_chol_mito']**(1.5),
                                                   lambda a : 0],
                    output_place_ids = ['p_H2O_mito'],  
                    production_speed_functions = [lambda a : k_t_ROS_metab * a['p_ROS_mito'] / a['p_chol_mito']**(1.5)])
                    

##Apoptosis
def add_apoptosis(hfpn):

    ## Places
    hfpn.add_place(initial_tokens=it_p_PERK_ER, place_id="p_PERK_ER", label="Conc of PERK in ER", continuous=True)

    hfpn.add_place(initial_tokens=it_p_BAX_mito, place_id="p_BAX_mito", label="Conc of BAX protein in ER", continuous=True)

    hfpn.add_place(initial_tokens=it_p_Bcl2_mito, place_id="p_Bcl2_mito", label="Conc of Bcl-2 regulator proteins in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_cytc_cyto, place_id="p_cytc_cyto", label="Conc of cytochrome C in mito", continuous=True)

    hfpn.add_place(initial_tokens=it_p_TRADD_cyto, place_id="p_TRADD_cyto", label="Conc of TRADD protein in cyto", continuous=True)

    hfpn.add_place(initial_tokens=it_p_mTORC1_LE, place_id="p_mTORC1_LE", label="Conc of mTORC1 in LE", continuous=True)

    hfpn.add_place(initial_tokens=it_p_cas3, place_id="p_cas3", label="Conc of caspase 3 in cyto", continuous=True)

    hfpn.add_place(initial_tokens=it_p_NFkB_cyto, place_id="p_NFkB_cyto", label="Conc of NF-kB proteins complex in cyto", continuous=True)

    # Places from other modules
    hfpn.add_place(initial_tokens=it_p_ROS_mito, place_id="p_ROS_mito", label="Conc of ROS in mito", continuous=True)
    hfpn.add_place(initial_tokens=it_p_tauP, place_id="p_tauP", label="Tau phosphorylated in cyto", continuous=True)
    hfpn.add_place(initial_tokens=it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_SINK, place_id="p_SINK", label="Sink for cholesterol export from PM", continuous=True)
    hfpn.add_place(initial_tokens=it_p_NPC1_LE, place_id="p_NPC1_LE", label="Mutated NPC1 in LE", continuous=True)
    hfpn.add_place(initial_tokens=it_p_Ca_LE, place_id="p_Ca_LE", label="Ca conc in LE", continuous=True)

    ## Transitions
    hfpn.add_transition(  transition_id = 't_ERstress_Ca_cyto_induc', 
                    label = 'ER stressed induced by Ca concentration in cyto', 
                    input_place_ids = ['p_ROS_mito', 'p_tauP','p_PERK_ER'],
                    firing_condition = fc_t_ERstress_Ca_cyto_induc,
                    consumption_speed_functions = [lambda a : 0,
                                                   lambda a : 0,
                                                   lambda a : 0],
                    output_place_ids = ['p_PERK_ER'],  
                    production_speed_functions = [lambda a : (k_t_ERstress_Ca_cyto_induc *(a['p_ROS_mito'] <= 23000)*9935.0758 * a['p_tauP'] + k_t_ERstress_Ca_cyto_induc *(a['p_ROS_mito'] > 23000)*a['p_ROS_mito'] * a['p_tauP'])*(a['p_PERK_ER'] < 10**8) + 10 * (a['p_PERK_ER'] >= 10**8)]) #12500 as average of healty ROS

    hfpn.add_transition(  transition_id = 't_BAX_actv', 
                    label = 'Activation of BAX', 
                    input_place_ids = ['p_PERK_ER', 'p_Bcl2_mito'],
                    firing_condition = fc_t_BAX_actv,
                    consumption_speed_functions = [lambda a : 10**5 * (a['p_PERK_ER'] < 10**8) + 10 * (a['p_PERK_ER'] >= 10**8), #Perk degradation rate
                                                   lambda a : k_t_Bcl2_actv * a['p_Bcl2_mito']**2 / it_p_Bcl2_mito**2], 
                    output_place_ids = ['p_BAX_mito'],  
                    production_speed_functions = [lambda a : k_t_BAX_actv * a['p_PERK_ER'] / (a['p_Bcl2_mito']+1) * 1/np.pi * (np.pi / 2 - np.arctan(a['p_Bcl2_mito'] - 4*10**5))])

    hfpn.add_transition(  transition_id = 't_Bcl2_actv', 
                    label = 'Activation of Bcl-2', 
                    input_place_ids = ['p_PERK_ER'],
                    firing_condition = fc_t_Bcl2_actv,
                    consumption_speed_functions = [lambda a : 0], 
                    output_place_ids = ['p_Bcl2_mito'],  
                    production_speed_functions = [lambda a : (a['p_PERK_ER'] < 10**7)*k_t_Bcl2_actv + (a['p_PERK_ER'] >= 10**7)*(k_t_Bcl2_actv * 10**3.5 / a['p_PERK_ER']**(1/2))])


    hfpn.add_transition(  transition_id = 't_mito_dysfunc', 
                    label = 'Dysfunction of mitochondrial Complex I', 
                    input_place_ids = ['p_BAX_mito', 'p_ROS_mito'],
                    firing_condition = fc_t_mito_dysfunc,
                    consumption_speed_functions = [lambda a : 33.33 * a['p_BAX_mito']**0.5,
                                                   lambda a : 0],
                    output_place_ids = ['p_cytc_cyto'],  
                    production_speed_functions = [lambda a : k_t_mito_dysfunc * a['p_BAX_mito']**0.5 + (a['p_ROS_mito'] > 23000) * k_t_mito_dysfunc_2 * 80000**0.5])

    hfpn.add_transition(  transition_id = 't_TRADD_actv', 
                    label = 'Activation of TRADD', 
                    input_place_ids = ['p_chol_LE'],
                    firing_condition = fc_t_TRADD_actv,
                    consumption_speed_functions = [lambda a : 0],
                    output_place_ids = ['p_TRADD_cyto'],  
                    production_speed_functions = [lambda a : k_t_TRADD_actv * a['p_chol_LE']**(0.2)])

    hfpn.add_transition(  transition_id = 't_Cas3_actv', 
                    label = 'Activtion of Caspase 3 apoptotic pathway', 
                    input_place_ids = ['p_TRADD_cyto', 'p_cytc_cyto'],
                    firing_condition = fc_t_Cas3_actv, 
                    consumption_speed_functions = [lambda a : k_t_Cas3_actv * a['p_TRADD_cyto'],
                                                   lambda a : k_t_Cas3_actv_2 * a['p_TRADD_cyto'] + (a['p_cytc_cyto'] >= 54228160) * 0.00977 * a['p_cytc_cyto']],
                    output_place_ids = ['p_cas3'],  
                    production_speed_functions = [lambda a : 8.1 * k_t_Cas3_actv * a['p_TRADD_cyto'] + (a['p_cytc_cyto'] >= 54228160) * 0.00977 * a['p_cytc_cyto']])

    hfpn.add_transition(  transition_id = 't_Cas3_degr', 
                    label = 'Degredation of Cas3', 
                    input_place_ids = ['p_cas3'],
                    firing_condition = fc_t_Cas3_degr,
                    consumption_speed_functions = [lambda a : k_t_Cas3_degr * (a['p_cas3']<=100000) + 2.7*10**6 * (a['p_cas3']>100000)], 
                    output_place_ids = ['p_SINK'],  
                    production_speed_functions = [lambda a : 0])

    hfpn.add_transition(  transition_id = 't_mTORC1_actv', 
                    label = 'Activation of mTORC1 cell proliferation pathway', 
                    input_place_ids = ['p_Ca_LE', 'p_NPC1_LE'],
                    firing_condition = fc_t_mTORC1_actv,
                    consumption_speed_functions = [lambda a : 0,
                                                   lambda a : 0],
                    output_place_ids = ['p_mTORC1_LE'],  
                    production_speed_functions = [lambda a : k_t_mTORC1_actv / a['p_NPC1_LE'] * ((a['p_Ca_LE'] > 220000)* (2.44*10**5)**(1/2) + (a['p_Ca_LE'] <= 220000)* a['p_Ca_LE']**(1/2))])

    hfpn.add_transition(  transition_id = 't_NFkB_actv', 
                    label = 'Actv of NF-kB inflam. signalling pathway', 
                    input_place_ids = ['p_mTORC1_LE', 'p_TRADD_cyto'],
                    firing_condition = fc_t_NFkB_actv,
                    consumption_speed_functions = [lambda a : k_t_NFkB_actv_2 * a['p_TRADD_cyto'] * a['p_mTORC1_LE'],
                                                   lambda a : k_t_NFkB_actv_3 * a['p_TRADD_cyto'] * a['p_mTORC1_LE']],
                    output_place_ids = ['p_NFkB_cyto'],  
                    production_speed_functions = [lambda a : k_t_NFkB_actv * a['p_TRADD_cyto'] * a['p_mTORC1_LE']])

    hfpn.add_transition(  transition_id = 't_NFkB_degr', 
                    label = 'Degredation of NFkB', 
                    input_place_ids = ['p_NFkB_cyto'],
                    firing_condition = fc_t_NFkB_degr, 
                    consumption_speed_functions = [lambda a : k_t_NFkB_degr],
                    output_place_ids = ['p_SINK'],  
                    production_speed_functions = [lambda a : 0])
