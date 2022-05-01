
#Matplotlib Graphing
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
style.use("ggplot")
import matplotlib.pyplot as plt

import numpy as np
#np.set_printoptions(threshold=100**3) ##jc/

import random


import pandas as pd

import tkinter as tk
from tkinter import ttk
from functools import partial
import math
import time

MAX_LABEL_CHARACTERS = 50

def check_label_length(label, limit = MAX_LABEL_CHARACTERS):
    if len(label) > limit:
            raise ValueError(f'Label \"{label}\" exceeds {MAX_LABEL_CHARACTERS} characters.')

#jc/ [class HFPN -> add_place() --> HFPN2.0 -> def PD_places()]
class Place:
    def __init__(self, initial_tokens, place_id, label, continuous = True):
        """
            Args:
                initial_tokens: number of tokens that the place is initialized with
                place_id (str): unique identifier of the place
                label (str): short description of the place
                continuous (bool): whether the place is continuous or discrete
        """
        self.place_id = place_id
        self.label = label
        self.initial_tokens = initial_tokens
        self.tokens = initial_tokens
        self.continuous = continuous

    def reset(self):
        self.tokens = self.initial_tokens

    def __str__(self):
        # TODO: update
        return f"{self.label} has {self.tokens} tokens"

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Place):
            return self.id == other.id
        return False


class ConsumptionSpeed: # The new input arc
    """The consumption speed describes at what rate tokens are consumed from its associated input place.
       Each input arc has an associated consumption speed, which generally depends on the concentration
       of each input place to the transition in question."""

    #BSL: downstream, consumption_place and input_places are identical, but they are passed differently in the functions contained within this class.
    def __init__(self, consumption_place, input_places, consumption_function):
        """
            Args:
                consumption_place: Place instance corresponding to the input place whose tokens are consumed
                input_places: list of Place instances corresponding to all the input places to the transition in question. 
                              These are necessary because the consumption speed generally depends on the concentration of
                              all the input places.
                consumption_function: lambda function which calculates the number of tokens to be consumed from consumption_place
        """
        self.consumption_place = consumption_place
        self.input_places = input_places
        self.consumption_function = consumption_function
        self.flag = 0
    
    #BSL:
    def __str__(self):
        return f"consumption_place:{self.consumption_place}, input_places:{self.input_places},consumption_function:{self.consumption_function}"
    
    def get_input_place_tokens(self):
        input_place_tokens_dictionary = {}
        for ip in self.input_places:
            input_place_tokens_dictionary[ip.place_id] = ip.tokens 
        return input_place_tokens_dictionary

    #BSL:
    def calculate_firing_tokens(self, time_step, stochastic_multiplier): #run each time fire() runs for a transition
        self.firing_tokens = self.consumption_function(self.get_input_place_tokens()) * time_step*stochastic_multiplier
        #if removing the no. tokens in firing_tokens will make place.tokens negative:
        if self.firing_tokens > self.consumption_place.tokens: #XXX we set -1 here, because this increases the WIP BROKEN
            self.flag =+ 1
    #BSL:
    def set_firing_tokens(self, token_value):
        self.firing_tokens = token_value
        
    def return_self_dot_firing_tokens(self):
        return self.firing_tokens
    
    def return_consumption_place_tokens(self):
        return self.consumption_place.tokens
    

    def perform_firing(self):
            self.consumption_place.tokens -= self.firing_tokens
            if self.consumption_place.tokens < 0:
                print(self.consumption_place.place_id)
                print(self.consumption_place.tokens)
                
    
   #BSL:         
    def return_consumed_fired_tokens(self): #maybe can improve by specifying storage intervals...
        return self.firing_tokens
    
    def reset_flag(self): #run every time after ive emptied a place of tokens, after detecting it being flagged
        self.flag = 0


class ProductionSpeed: # The new output arc
    """The production speed describes at what rate tokens are produced to its associated output place.
       Each output arc has an associated production speed, which generally depends on the concentration
       of each input place to the transition in question."""

    def __init__(self, production_place, input_places, production_function):
        """
            Args:
                production_place: Place instance corresponding to the output place where tokens are produced.
                input_places (list): list of Place instances corresponding to all the input places to the transition in question. 
                                     These are necessary because the production speed generally depends on the concentration of
                                     all the input places.
                production_function: lambda function which calculates the number of tokens to be produced in production_place
        """
        self.production_place = production_place
        self.input_places = input_places
        self.production_function = production_function

    #BSL:
    def __str__(self):
        return f"production_place:{self.production_place}, input_places:{self.input_places},production_function:{self.production_function}"
    
    def get_input_place_tokens(self):
        input_place_tokens_dictionary = {}
        for ip in self.input_places:
            input_place_tokens_dictionary[ip.place_id] = ip.tokens 
        #print(input_place_tokens_dictionary)
        return input_place_tokens_dictionary


    def calculate_firing_tokens(self, time_step, stochastic_multiplier):
        self.firing_tokens = self.production_function(self.get_input_place_tokens()) * time_step * stochastic_multiplier
        
        #print(self.firing_tokens, "production firing_tokens")
    #BSL:
    def set_firing_tokens(self, token_value):
        self.firing_tokens = token_value

    def perform_firing(self):
        self.production_place.tokens += self.firing_tokens 
   #BSL: 
    def return_produced_fired_tokens(self):
        return self.firing_tokens

class ContinuousTransition:
    """A continuous transition contains (i) a firing condition, usually expressed in terms of the input concentrations,
       (ii) the collection of all consumption speeds, and (iii) the collection of all production speeds."""

    def __init__(self, 
                 transition_id, 
                 label, 
                 input_place_ids, 
                 firing_condition, 
                 consumption_speeds, 
                 production_speeds, 
                 stochastic_parameters,
                 collect_rate_analytics, 
                 consumption_coefficients, 
                 output_place_ids,
                 production_coefficients):
        """
            Args:
                transition_id (str): unique identifier of a transition
                label (str): short description of the transition
                firing_condition: lambda function which describes the criteria for the transition to fire
                consumption_speeds (list): each element is the consumption speed from an input place
                production_speeds (list): each element is the production speed to an output place
        """
        self.transition_id = transition_id
        self.label = label
        self.input_place_ids=input_place_ids
        self.firing_condition = firing_condition
        self.consumption_speeds = consumption_speeds
        self.production_speeds = production_speeds
        self.stochastic_parameters = stochastic_parameters
        self.collect_rate_analytics = collect_rate_analytics
        self.firings = 0
        self.list_of_consumed_tokens = []
        self.list_of_produced_tokens = []
        self.consumption_coefficients = consumption_coefficients #BSL
        self.production_coefficients = production_coefficients
        self.counter_thing = 0
        self.output_place_ids = output_place_ids
        self.DiscreteFlag = "no"
        
    #BSL:
    def __str__(self):
        return f"transition_id:{self.transition_id}, label:{self.label},firing_condition:{self.firing_condition},consumption_speeds:{self.consumption_speeds},production_speeds:{self.production_speeds},firings:{self.firings}"
    
    

    def get_input_place_tokens(self):
        input_place_tokens_dictionary = {}
        for cs in self.consumption_speeds:
            input_place_tokens_dictionary[cs.consumption_place.place_id] = cs.consumption_place.tokens 
        return input_place_tokens_dictionary

    def reset(self):
        """Reset firing count for each transiton to 0 after each run."""
        self.firings = 0
    def reset_list_of_produced_tokens(self):
        self.list_of_produced_tokens = []
        
    def reset_list_of_consumed_tokens(self):
        self.list_of_consumed_tokens = []

    def fire(self, time_step):
        """The function first checks whether the firing condition is satisfied. If it is, it first
           calculates the amount tokens that will be transferred/produced from/to each place. Finally,
           this transfer is actually performed.
        """
        input_place_tokens = self.get_input_place_tokens()
            #= dict: key = input_place_id, val = tokens
            
        # Check if the firing condition of the transition is satisfied
        if self.firing_condition(input_place_tokens) == True:
            # Calculate number of tokens that will be consumed or produced from firing
            #calculated tokens list
            #randomizer
            randomized_value=random.gauss(1,self.stochastic_parameters[0])
            for cs in self.consumption_speeds:
                cs.calculate_firing_tokens(time_step, randomized_value)
            
            #Check if flags are made
            check_flag = np.array([]) 

            for cs in self.consumption_speeds:
                check_flag = np.append(check_flag, cs.flag) #delete check_flag and u will see negative values appearing
                #jc/one by one, append the flag boolean onto the check_flag array for each input place to the transition
                    #(eg if 3 places in {} input_place_tokens -> 3 (0/1) flags,
                    #that are appended 1 by 1, to produce a 1x3 check_flag array)

#jc/fire transition: if flags ------------------------------------------------------------------
            if any(check_flag): #skips this block if no flags

                calculated_tokens_list = np.array([])
                # print("")#DEBUGGING
                # print(self.transition_id) #DEBUGGING   
                # print("")#DEBUGGING
                for cs in self.consumption_speeds: #consumption_speeds = list of instances of ConsumptionSeed(), for each place listed in input_place_ids for the transition im running in fire()
                    # print("")#DEBUGGING
                    # print(cs.firing_tokens, "pre-correction FIRING TOKENS")#DEBUGGING
                    # print("")#DEBUGGING
                    if cs.flag == 0: #cs = instance of ConsumptionSpeed()
                        calculated_tokens_list = np.append(calculated_tokens_list, np.nan)
                    if cs.flag >0: #if flagged, then Empty place contents and fire tokens = limiting factor
                        calculated_tokens_list = np.append(calculated_tokens_list,cs.return_consumption_place_tokens())
                        # print("input place flagged: " + str(cs.consumption_place.place_id) + " (" + str(cs.consumption_place.tokens)+" tokens)")
                #divide calculated tokens by consumption coefficients to identify limiting factors from a list of multiple candidate limiting factors.     
                deciding_value_list = [i/j for i,j in zip(calculated_tokens_list,self.consumption_coefficients)]#eg [nan, 10, nan] or [13, 10, nan] #dividing the consumption_place_tokens by the consumption coefficients should tell you which element in the list to prioritise.
                #select then scale
                index_min = np.nanargmin(deciding_value_list) #select index with the lowest consumption_tokens, and thus is the minimum value, this is the limiting factor. All the tokens for this prioritised list element should be completely emptied. index_min is an integer, typically 0, 1 or 2.
                #Need to produce the list of tokens to be consumed:
                current_place_tokens = np.array([])
                for cs in self.consumption_speeds:
                    current_place_tokens = np.append(current_place_tokens, cs.return_consumption_place_tokens())
                list_len = len(current_place_tokens)
                consuming_tokens = np.zeros(list_len)
                tokens_of_prioritised_element = (current_place_tokens[index_min]) -1 # minus 1 is crucial to avoid divide by zero errors,so this empties the prioritised element to 1 token left, instead of 0 tokens left. #BROKEN
                standardised_tokens_of_prioritised_element = (tokens_of_prioritised_element/self.consumption_coefficients[index_min])
                for index,coefficient in enumerate(self.consumption_coefficients): #consumption_coefficients list should match to length of consuming_tokens list
                    consuming_tokens[index] = standardised_tokens_of_prioritised_element*coefficient
                
                #Now we need to set the tokens to fire for consumption_speeds
                for cs,token_value in zip(self.consumption_speeds, consuming_tokens):
                    token_value = token_value 
                    cs.set_firing_tokens(token_value) #
                
                #calculate produced tokens from production_coefficients
                list_len2 = len(self.production_coefficients)
                producing_tokens = np.zeros(list_len2)
                for index,coefficient in enumerate(self.production_coefficients): 
                    producing_tokens[index] = standardised_tokens_of_prioritised_element*coefficient                
                    
                #Now we need to set the tokens to be produced for production_speeds:
                for ps, token_value_prod in zip(self.production_speeds, producing_tokens):
                    ps.set_firing_tokens(token_value_prod)
                
                #Perform Firing
                for cs in self.consumption_speeds:
                    # print(cs.get_input_place_tokens(),"input place tokens Before Firing post correction") #DEBUGGING
                    cs.perform_firing()
                    if self.collect_rate_analytics[0]=="yes":
                        self.list_of_consumed_tokens.append(cs.return_consumed_fired_tokens())
    
                for ps in self.production_speeds:
                    # print(ps.production_place.tokens,"production_place tokens before firing") #DEBUGGING
                    ps.perform_firing()
                    if self.collect_rate_analytics[1] == "yes":
                        self.list_of_produced_tokens.append(ps.return_produced_fired_tokens())

                for cs in self.consumption_speeds: #important to reset flag 
                    cs.reset_flag()
                #debugging chunk start
                # print(self.transition_id, consuming_tokens, "consumingTokens")
                # print(self.consumption_coefficients, "consumption_coefficients")
                # print(self.production_coefficients, "production_coefficients")
                # print(self.consumption_speeds[0].get_input_place_tokens(), "after firing post correction")
                # for index,ps in enumerate(self.production_speeds):
                    # print(ps.production_place.tokens, "production Place tokens after firing", index+1)
                    # print(ps.return_produced_fired_tokens(), "produced fired tokens", index+1)
                    # print("")#DEBUGGING
                    # print(randomized_value, "Randomized Value")
                self.counter_thing +=1
                #print(self.counter_thing, "flag counter for", self.transition_id)
                #debugging end
                
#jc/fire transition: if NO flags ------------------------------------------------------------------
            else: #if check flag is a list of pure zeros, then this block of code should get executed.
                
                for ps in self.production_speeds:
                    ps.calculate_firing_tokens(time_step, randomized_value)
     
                # store least non-zero transfer count
                token_transfers = [s.firing_tokens for s in self.consumption_speeds if s.firing_tokens != 0] + \
                                  [s.firing_tokens for s in self.production_speeds if s.firing_tokens != 0]
                                
                if len(token_transfers) == 0: 
                    token_transfers.append(0)
                
                # Execute the actual firing by looping through all ConsumptionSpeed and ProductionSpeed instances            
                for cs in self.consumption_speeds: #jc/cs = instance of ConsumptionSpeed()
                    #jc/ for 1 input place for 1 transition:
                    cs.perform_firing() #jc//TAKE AWAY TOKENS FROM INPUT PLACES
                    if self.collect_rate_analytics[0]=="yes": #self = instance of ContinuousTransition (the transition currently being fired)
                        # print("========entered fire()") #jc/CS_CLOCK
                        # print(self.transition_id)
                        # print(cs.input_places[0].place_id)
                        # print(cs.consumption_place.place_id)
                        # print("==============================")
                        self.list_of_consumed_tokens.append(cs.return_consumed_fired_tokens())
    
                for ps in self.production_speeds:
                    ps.perform_firing() #jc//ADD TOKENS TO OUTPUT PLACES - CHANGES NO. TOKENS
                #BSL:
                #store produced tokens PER TRANSITION for later analysis #this is wrong if more than one thing produced I think
                #self.list_of_consumed_tokens.append(return_consumed_fired_tokens())
                    if self.collect_rate_analytics[1] == "yes":
                        self.list_of_produced_tokens.append(ps.return_produced_fired_tokens())
            
            # Increment number of firings by 1
            self.firings += 1
        #this is so that we know the rate is zero, if the transition doesn't fire in rate analytics. We need to change this to a numpy array. #need to also do the same for list_of_consumed_tokens when we get to it.
        
#jc/transition does NOT meet firing requirements -----------------------------------        
        #If there are two output places, then 2 zeros are appended into the list.
        if self.firing_condition(input_place_tokens) == False:

            for cs in self.consumption_speeds:
                if self.collect_rate_analytics[0]=="yes":
                    self.list_of_consumed_tokens.append(0)
            for ps in self.production_speeds:
                if self.collect_rate_analytics[1] =="yes":
                    self.list_of_produced_tokens.append(0)
            

            
class DiscreteTransition(ContinuousTransition):

    def __init__(self, 
                 transition_id, 
                 label,input_place_ids, 
                 firing_condition, 
                 consumption_speeds, 
                 production_speeds,
                 stochastic_parameters,
                 collect_rate_analytics, 
                 delay, 
                 consumption_coefficients, 
                 output_place_ids, 
                 production_coefficients):
        """In addition to the arguments specified in the super class ContinuousTransition, a delay function must 
           be specified for a discrete transition.

            Args:
                delay (int): number of time-steps after which transition is fired if firing_condition still holds true.
        """

        # Initialize everything from the super class
        super(DiscreteTransition, self).__init__(transition_id, 
                                                 label,
                                                 input_place_ids, 
                                                 firing_condition, 
                                                 consumption_speeds,
                                                 production_speeds,
                                                 stochastic_parameters,
                                                 collect_rate_analytics, 
                                                 consumption_coefficients,
                                                 output_place_ids,  
                                                 production_coefficients)
        self.delay = delay
        self.delay_original = delay
        self.delay_counter = 0 # Counter for consecutive no. of steps where firing condition still holds true
        self.delay_list = [] # create list for analysis later
        self.DiscreteFlag = "yes"
        
        
        input_place_tokens = self.get_input_place_tokens()
        self.delay = self.calculate_lambda_f_delay(self.delay_original, input_place_tokens)

    def calculate_lambda_f_delay(self, delay_original, input_place_tokens):
        """Purpose is to check if delay is a lambda function or a float/integer, only used in initialisation to prevent error"""
        truth = isinstance(delay_original, int)
        truth2 = isinstance(delay_original, float)
        
        if int(truth) + int(truth2) == 0:
            delay = self.delay_original(input_place_tokens)
            self.delay_original
            return delay
            
        else:
            return delay_original
            
        


    def fire(self, time_step):
        """Check if the firing condition is satisfied during the delay."""
        input_place_tokens = self.get_input_place_tokens()
        if self.firing_condition(input_place_tokens) == True:
            # if self.transition_id == "t_MDV_Generation_basal": #debugging
            #     print(self.transition_id)
            #     print(self.delay)
            #     print(self.delay_counter)
            #     print(int(self.delay/time_step))

            if self.delay_counter >= int(self.delay/time_step):
                # fire with a time_step of 1, as discrete transition tokens should be transferred in their entirety
                super().fire(1) #BSL: Means you are calling fire() from the Parent Class ContinuousTransition, from parent class fire, you should have appended consumed/produced tokens for rate analytics.
                self.delay_list.append(self.delay_counter) #append to this list for later analysis
                self.delay_counter = 0
                
                if isinstance(self.delay_original, (int,float)):
                    # print("---start------------------")
                    # print(self.transition_id)
                    # print(self.delay)
                    # print(self.delay_counter)
                    # print("--------------------")
                    # print(self.stochastic_parameters)
                    # #print(self.stochastic_parameters[1])
                    # print("---end------------------")
                    self.delay = random.gauss(self.delay_original,self.delay_original*self.stochastic_parameters[1])
                else:
                    value = self.delay_original(input_place_tokens)
                    self.delay = random.gauss(value, value*self.stochastic_parameters[1])
                    
            else:
                self.delay_counter += 1
                for ps in self.production_speeds:
                    if self.collect_rate_analytics[1]=="yes":
                        self.list_of_produced_tokens.append(0)   
                        
                for cs in self.consumption_speeds:
                    if self.collect_rate_analytics[0]=="yes":
                        self.list_of_consumed_tokens.append(0)
        else:
            # Reset firing condition
            self.delay_counter = 0
            for ps in self.production_speeds:
                if self.collect_rate_analytics[1]=="yes":
                    self.list_of_produced_tokens.append(0)  
            for cs in self.consumption_speeds:
                if self.collect_rate_analytics[0]=="yes":
                    self.list_of_consumed_tokens.append(0)
            
#jc/ class HFPN [-->HFPN2.0]

class HFPN:
    """Hybrid functional Petri net (HPFN) with the option to specify a time-step in seconds."""

    def __init__(self, printout=False): 
        """
            Args:
                printout (bool): Whether to print out number of tokens in each place for each time-step
        """
        self.places = {} #jc/instance variables aka attributes - set as empty dict
        self.transitions = {}
        self.printout = printout #jc/ie = attri val, default = False
        self.counter = 0
        
#jc/ set_place_tokens() [<- add_place()]

    def set_place_tokens(self, place_id, value):
        self.places[place_id].tokens = value
        #1. attri self.places = dict, containing key-val of every place
            #val: instances of class Place()
        #2. self.places["p_chol_PM"] = instance of class Place()
            #contains the attri initial_tokens, label, etc, pertaining to "p_chol_PM"
        #3. self.places["p_chol_PM].tokens = attri (self.tokens)
            #ie tokens = attri of class instance self.places["p_chol_PM]
            #tokens attri stores val 'value'
            
    def remove_Ca_clock(self):
        # self.transitions['t_Ca_imp'].input_place_ids = []
            #self.transitions['t_Ca_imp'] = instance
        self.transitions['t_Ca_imp'].firing_condition = lambda a: True
        self.transitions['t_Ca_imp'].reaction_speed_function = lambda a : 1.44*1e8*1/3
        # self.transitions['t_Ca_imp'].consumption_coefficients = []
        
        # self.transitions['t_NCX_PMCA'].input_place_ids = ['p_Ca_cyto']
        self.transitions['t_NCX_PMCA'].firing_condition = lambda a: True
        self.transitions['t_NCX_PMCA'].reaction_speed_function = lambda a: (10*a['p_Ca_cyto'])*(a['p_Ca_cyto'])*1/3
        # self.transitions['t_NCX_PMCA'].consumption_coefficients = [1]
        
        # self.transitions['t_NaK_ATPase'].input_place_ids = ['p_ATP']
        self.transitions['t_NaK_ATPase'].firing_condition = lambda a: True
        self.transitions['t_NaK_ATPase'].rate_constant =  0.7*1/3
        # self.transitions['t_NaK_ATPase'].consumption_coefficients = [1]
        
        # self.transitions['t_RyR_IP3R'].input_place_ids = ['p_Ca_ER']
        self.transitions['t_RyR_IP3R'].firing_condition = lambda a: True
        self.transitions['t_RyR_IP3R'].reaction_speed_function = lambda a: (10*a['p_Ca_cyto'])*(a['p_Ca_cyto'])*1/3
        # self.transitions['t_RyR_IP3R'].consumption_coefficients =[1]
            
#jc/ set_1st_stochastic_parameter() [<- add_transition()]

    def set_1st_stochastic_parameter(self, decimal, transition_id):
        self.transitions[transition_id].stochastic_parameters[0] = decimal
        #1. attri self.transitions = DICT, containing key-val of every transition
            #val: instances of class DiscreteTransition()/ContinuousTransition()
        #2. self.transitions["t_CYP7B1_metab"] = INSTANCE of class DiscreteTransition()
            #contains the attri transition_id, label, etc, pertaining to "t_CYP7B1_metab"
        #3. #self.transitions["t_CYP7B1_metab"].stochastic_parameters 
                #= ATTRI of instance (self.stochastic_parameters)
                #= list, eg [1, 0, 0]
            #stochastic_parameters = attri of DiscreteTransition() class
            #inputs that are put into .stochastic_parameters = list, eg [1, 0, 0]
        #4. self.transitions[transition_id].stochastic_parameters[0]
                #= replaces the value of the 1st val in list 'stochastic_parameters' to 
                #the val of variable 'decimal'
            #(if .stochastic_parameters hadnt existed as an attri of instance 
            #transitions["t_CYP7B1_metab"], one would be newly created and set to
            #val 'decimal')
        #if i want to recall the value of decimal:
            #instance = HFPN()
            #instance.transitions["t_CYP7B1_metab"].stochastic_parameters[0]
            
    def set_2nd_stochastic_parameter(self, decimal, transition_id):
        "Allows the stochastic Parameter to be reset later, after the transition was loaded"
        self.transitions[transition_id].stochastic_parameters[1]=decimal
        

    def set_consumption_collect_decision(self, integer, transition_id):

        if integer == 0:
            self.transitions[transition_id].collect_rate_analytics[0] = "no"
            #print("0 activated")
            #print(self.transitions[transition_id].collect_rate_analytics, "IN PETRI NET")            
        else:
            self.transitions[transition_id].collect_rate_analytics[0] = "yes"
            #print("1 activated")
            #print(self.transitions[transition_id].collect_rate_analytics, "IN PETRI NET")
            
        
    def set_production_collect_decision(self, integer, transition_id):
        if integer == 0:
            self.transitions[transition_id].collect_rate_analytics[1] = "no"
        if integer ==1:
            self.transitions[transition_id].collect_rate_analytics[1] = "yes"
   
        
    def set_time_step(self, time_step):
        "time_step: time increment size (seconds)"
        self.time_step = time_step
            #BSL:
    def __str__(self):
        return f"places:{self.places}, transitions:{self.transitions},time_step:{self.time_step}"

    def set_initial_tokens_for(self, place_id, initial_tokens):
        self.places[place_id].initial_tokens = initial_tokens
        self.places[place_id].tokens = initial_tokens

#jc/ add_place() [--> HFPN2.0 -> def PD_places()]
    def add_place(self, initial_tokens, place_id, label, continuous = True):
        """Add a place to the Petri net.

            Args:
                initial_tokens (int): initial number of tokens in this place
                place_id (str): unique identifier of this place
                label (str): description of the place
                continuous (bool): whether the place is continuous or discrete
        """ #???jc/^what does contin/discrete mean
        check_label_length(label) # Check whether the label is too long

        # Check whether the place_id contains a space
        if ' ' in place_id:
            raise ValueError(f"Place_id should not contain any spaces {place_id}. Did you reverse place_id and label?")

        if place_id in self.places.keys():
            print(f"Warning: Place {place_id} already exists.")
            #rewrite it.. does this work?
            place = Place(initial_tokens, place_id, label, continuous)
            self.places[place_id] = place
            print(f"Rewrote: Place {place_id}.")
        else:
            place = Place(initial_tokens, place_id, label, continuous)
            #jc/create new instance of Place() class
                #initial_tokens, place_id, etc become 'place' attributes 'place.attri'
            self.places[place_id] = place
            #jc/ 1. self.places = empty dict {}
                #2. for EACH new add_place() called: create new key-value pair in places{}
                    #key: 'place_id' (= "str")
                    #value: 'place' instance of class Place(), containing info about the place
                #3. attri self.places = refers computer to the dictionary object
                    #dict obj = contains key-val pair of every place

#jc/ add_transition() [--> HFPN2.0 -> def ()]
    def add_transition( self,
                        transition_id,
                        label,
                        input_place_ids,
                        firing_condition,
                        consumption_speed_functions,
                        output_place_ids,
                        production_speed_functions,
                        stochastic_parameters,
                        collect_rate_analytics,
                        consumption_coefficients,
                        production_coefficients,
                        delay = -1):
        """Adds a transition to the hybrid functional Petri net.

            Args: 
                transition_id (str): unique identifier of the transition
                label (str): short description of the transition
                input_place_ids (list): list with place_id for each input place
                firing_condition: lambda function which returns a Boolean
                consumption_speed_functions (list): list of lambda functions, each of which returns a float or int
                output_place_ids (list): list with place_id for each output place
                production_speed_functions (list): list of lambda functions, each of which returns a float or int
                delay (float): the delay, in seconds, associated with a discrete transition. If not specified, a continuous transition is assumed.
        """
        # Check if function inputs are functions
        for function in consumption_speed_functions+production_speed_functions+[firing_condition]:
            if not callable(function): #BSL: lambda functions are functions, and functions are callable
                raise TypeError(f"Transition {transition_id} contains a non-function consumption speed function/production speed function or firing condition: {function}")
        # Check if input/output places exist
        for place_id in input_place_ids+output_place_ids: 
            if place_id not in self.places.keys():
                raise ValueError(f"Place {place_id} in transition {transition_id} has not been defined")
        # Check that the parameters of the function have been passed properly
        if len(input_place_ids) != len(consumption_speed_functions):
            raise ValueError(f"Unequal numbers of input places and input arc weights in transition {label}")
        if len(output_place_ids) != len(production_speed_functions):
            raise ValueError(f"Unequal numbers of output places and output arc weights in transition {label}")

        check_label_length(label) # Check whether the label is too long

        # Check whether the transition_id contains a space
        if ' ' in transition_id:
            raise ValueError(f"Transition_id should not contain any spaces {transition_id}. Did you reverse transition_id and label?")

        if transition_id in self.transitions.keys():
            print(f"Warning: Transition {transition_id} already exists. Old {transition_id} values have been replaced")
            #Rewrite
            
        # else:

        # Translate input_places from strings to Place instances
        consumption_speeds = []
        for ipid, csf in zip(input_place_ids, consumption_speed_functions): #so each object in consumption_speeds is an input_place_id key with the corresponding consumption_speed_function
            consumption_speeds.append(ConsumptionSpeed(self.places[ipid], self.places_from_keys(input_place_ids), csf)) #first two are identical, because you need to pass the consumption place argument AND the input place argument
                        
        production_speeds = []
        for opid, psf in zip(output_place_ids, production_speed_functions):
            production_speeds.append(ProductionSpeed(self.places[opid], self.places_from_keys(input_place_ids), psf))

        #jc/ *main segment*
        if delay == -1: # user wants continous transition as delay is not specified
            transition = ContinuousTransition(transition_id, label, input_place_ids,firing_condition, consumption_speeds, production_speeds, stochastic_parameters,collect_rate_analytics, consumption_coefficients,output_place_ids, production_coefficients)
                #transition = instance of ContTran()
            # Check if continuous transition  is linked to discrete output place, which is not allowed.
            output_places = self.places_from_keys(output_place_ids)
            for op in output_places:
                if op.continuous == False:
                    raise ValueError(f"A continuous transition ({transition_id}) cannot be linked to a discrete output place ({op.place_id}).")             
        else:         
            transition = DiscreteTransition(transition_id, label,input_place_ids, firing_condition, consumption_speeds, production_speeds, stochastic_parameters,collect_rate_analytics, delay, consumption_coefficients,output_place_ids, production_coefficients) 
            #jc/create new instance of DiscreteTransition() class
                #transition_id, label, etc get fed in and become 'transition' attributes 'transition.attri'
        self.transitions[transition_id] = transition
        #jc/ 1. self.transitions = empty dict {}
            #2. for EACH new add_transition() called: create new key-value pair in places{}
                #key: 'transition_id' (= "str")
                #value: 'transition' instance of class DiscreteTransition(), containing info about the transition
            #3. attri self.transitions = refers computer to the dictionary object self.transitions
                #dict obj = contains key-val pair of every transition  
                           
#jc/ add_transition_with_speed_function
    def add_transition_with_speed_function(self,
                        transition_id,
                        label,
                        input_place_ids,
                        firing_condition,
                        reaction_speed_function, ##
                        consumption_coefficients,
                        output_place_ids,
                        production_coefficients,
                        stochastic_parameters,
                        collect_rate_analytics,
                        delay = -1):
        
        """Add transition to HFPN wherein all consumption/production speeds are defined as proportional to a given reaction_speed_function.

            Args: 
                transition_id (str): unique identifier of the transition
                label (str): short description of the transition
                input_place_ids (list): list with place_id for each input place
                firing_condition: lambda function which returns a Boolean
                reaction_speed_function: lambda function which returns the number of tokens to transfer with dependence on input_places
                consumption_coefficients (list): list of numbers relating the consumption speed of the input places
                output_place_ids (list): list with place_id for each output place
                production_coefficients (list): list of numbers relating the production speed of the output places
                delay (float): the delay, in seconds, associated with a discrete transition. If not specified, a continuous transition is assumed.
        """
        # Nested lambda function, see https://www.geeksforgeeks.org/nested-lambda-function-in-python/
        for function in [reaction_speed_function, firing_condition]:
            if not callable(function):
                raise TypeError(f"Transition {transition_id} contains a non-function reaction speed function or firing condition: {function}")
        function = lambda f, n : lambda a : f(a) * n  #BSL: you are multiply by n because n is the stoichiometry of the reaction.

        consumption_speed_functions = [function(reaction_speed_function, cc) for cc in consumption_coefficients]
        production_speed_functions  = [function(reaction_speed_function, pc) for pc in production_coefficients] #BSL: List whose length corresponds to the number of production coefficients. 

        self.add_transition(
            transition_id = transition_id,
            label = label,
            input_place_ids = input_place_ids,
            firing_condition = firing_condition,
            consumption_speed_functions = consumption_speed_functions,
            output_place_ids = output_place_ids,
            production_speed_functions = production_speed_functions,
            stochastic_parameters=stochastic_parameters,
            consumption_coefficients = consumption_coefficients,
            production_coefficients = production_coefficients,
            collect_rate_analytics=collect_rate_analytics,
            delay = delay)
        

    def add_transition_with_mass_action( self,
                        transition_id,
                        label,
                        rate_constant,
                        input_place_ids,
                        firing_condition,
                        consumption_coefficients, 
                        output_place_ids,
                        production_coefficients,
                        stochastic_parameters,
                        collect_rate_analytics,
                        delay = -1):

        """Adds a transition to the HFPN where all firing rates are defined based on mass action.

            Args: 
                transition_id (str): unique identifier of the transition
                label (str): short description of the transition
                rate_constant (float): number multiplied by mass action rate equation to produce overall reaction speed 
                input_place_ids (list): list with place_id for each input place
                firing_condition: lambda function which returns a Boolean
                consumption_coefficients (list): list of numbers relating the consumption speed of the input places
                output_place_ids (list): list with place_id for each output place
                production_coefficients (list): list of numbers relating the production speed of the output places 
                delay (float): the delay, in seconds, associated with a discrete transition. If not specified, a continuous transition is assumed.
        """

        mass_action_function = lambda a : rate_constant * np.prod([a[place_id] ** ratio for place_id, ratio in zip(input_place_ids, consumption_coefficients)])

        self.add_transition_with_speed_function(
            transition_id = transition_id,
            label = label,
            input_place_ids = input_place_ids,
            firing_condition = firing_condition,
            reaction_speed_function = mass_action_function,
            consumption_coefficients = consumption_coefficients,
            output_place_ids = output_place_ids,
            production_coefficients = production_coefficients,
            stochastic_parameters = stochastic_parameters,
            collect_rate_analytics = collect_rate_analytics,
            delay = delay)


    def add_transition_with_michaelis_menten(    self,
                                                transition_id,
                                                label,
                                                Km,
                                                vmax,
                                                input_place_ids,
                                                substrate_id,
                                                consumption_coefficients,
                                                output_place_ids,
                                                production_coefficients,
                                                stochastic_parameters,
                                                collect_rate_analytics,
                                                vmax_scaling_function = (lambda a : 1)):
        """Adds a transition to the HFPN where firing rates are defined based on Michaelis Menten.

            Args: 
                transition_id (str): unique identifier of the transition
                label (str): short description of the transition
                Km (float): Michaelis Menten Km, substrate concentration for speed to be vmax/2
                vmax (float): Michaelis Menten vmax, maximum speed transition can achieve
                input_place_ids (list): list with place_id for each input place
                substrate_id (str): id of substrate of reaction
                consumption_coefficients (list): list of numbers relating the consumption speed of the input places
                output_place_ids (list): list with place_id for each output place
                production_coefficients (list): list of numbers relating the production speed of the output places
                vmax_scaling_function: lambda function specifying how promoters/inhibitors affect vmax
        """

        speed_function = lambda a : vmax * a[substrate_id] / (Km + a[substrate_id])
        function = lambda f, g : lambda a : f(a) * g(a) # compound the two funtions
        scaled_speed_function = function(speed_function, vmax_scaling_function)


        self.add_transition_with_speed_function(
            transition_id = transition_id,
            label = label,
            input_place_ids = input_place_ids,
            firing_condition = lambda a : True,
            reaction_speed_function = scaled_speed_function,
            consumption_coefficients = consumption_coefficients,
            output_place_ids = output_place_ids,
            production_coefficients = production_coefficients,
            stochastic_parameters = stochastic_parameters,
            collect_rate_analytics= collect_rate_analytics)


    
    def places_from_keys(self, keys): 
        #print(keys, "Keys") #debug
        #for key in keys:
            #print(self.places[key], "self.places[key]") #each component of the dictionary is a place, so its accessing __str__ of Place
        return [self.places[key] for key in keys]

    def find_places_transitions(self, string, case_sensitive = True, search_places=True, search_transitions=True):
        if search_places == True:
            print(f"List of all place id's containing {string}:")
            for pid in self.places.keys():
                if case_sensitive == True:
                    if string in pid:
                        print(f"{pid}")
                else: 
                    if string.lower() in pid.lower():
                        print(f"{pid}")

        if search_transitions== True:
            print(f"\nList of all transitions id's containing {string}:")
            for tid in self.transitions.keys():
                if case_sensitive == True:
                    if string in tid:
                        print(f"{tid}")
                else: 
                    if string.lower() in tid.lower():
                        print(f"{tid}")

    def find_places_transitions_labels(self, string, case_sensitive = True, search_places=True, search_transitions=True):
        if search_places == True:
            print(f"List of all place labels containing {string}:")
            for place in self.places.values():
                if case_sensitive == True:
                    if string in place.label:
                        print(f"{place.label}")
                else:
                    if string.lower() in place.label.lower():
                        print(f"{place.label}")

        if search_transitions== True:
            print(f"\nList of all transitions labels containing {string}:")
            for transition in self.transitions.values():
                if case_sensitive == True:
                    if string in transition.label:
                        print(f"{transition.label}")
                else:
                    if string.lower() in transition.label.lower():
                        print(f"{transition.label}")


    def run_single_step(self):
        """For each time-step, generate random order of all transitions and execute sequentially.

            Returns:
                tokens (list): list with number of tokens in each place
                firings (list): list with number of firings for each transitions
        """

        ordered_transitions = list(self.transitions.values())
            #jc/ .values() = inbuilt dict method that outputs just the values of
                #a dict (ie not the keys)
# =============================================================================
#         print("ordered_transitions-----------------------------------")
#         print(ordered_transitions) #jc/--> list of 51 instances from 2 dif classes:
#                                             45 x "stochastic_hfpn.ContinuousTransition object"
#                                             6 x "stochastic_hfpn.DiscreteTransition object"
#         print(len(ordered_transitions))
#         print(type(self.transitions.values())) #jc/--> dict_values
#         print(self.transitions.values()) #jc/--> collection of object instances of Continuous/DiscreteTransition()
# =============================================================================
        random_order_transitions = random.sample(ordered_transitions, len(ordered_transitions))
# =============================================================================
#         print("random_order_transitions--------------------------------")
#         print(random_order_transitions) #jc/--> list like ordered_transitions, but with the 51
#                                             items in randomised order
#         print(random_order_transitions[0]) #jc/--> = instance of class ContinuousTransition().
#                 Attri:
#                     transition_id:t_MAM, 
#                     label:Ca transport from ER to mitochondria
#                     firing_condition:<function <lambda> at 0x00000209A1CD8D30>
#                     consumption_speeds:[<stochastic_hfpn.ConsumptionSpeed object at 0x0000020B4E3FF370>, <stochastic_hfpn.ConsumptionSpeed object at 0x0000020B4E3FF3D0>]
#                     production_speeds:[<stochastic_hfpn.ProductionSpeed object at 0x0000020B4E3FF430>]
#                     firings:3204
#         print(random_order_transitions[0].label) #jc/--> [eg] Ca transport from ER to mitochondria
#         print(random_order_transitions[0].firings) #jc/--> no. betw ~0 & 4000
# =============================================================================
        
        for t in random_order_transitions:
            #jc/ each iteration: t becomes ONE of the 51 instances of class 
                #DiscreteTransition()/ContinuousTransition() stored within random_order_transitions
                #ie. t now has DiscreteTransition()/ContinuousTransition() attri 
                #t.transition_id, t.firing, etc and method t.fire()
            #jc/ every transition fires once per timestep (ie each timestep contains 51 firings)
            t.fire(self.time_step)

# =============================================================================
#             print("transition: " + str(t.transition_id) + "------------------------------------")
#             print("self.places.values())[10]:")
#             print(list(self.places.values())[10].tokens)
#                   #jc/--> it's always the same transition in random_order_transitions
#                       (t_24OHchol_exocyto) that changes the .tokens number of 
#                       the place in [10] (24OHchol)
# =============================================================================
            
        # Store tokens weights for each place at specific time step 
        tokens = [place.tokens for place in self.places.values()]
            #jc/ 1. self.places = {} dict:
            		#Key: initial_tokens
                    #Val: instance of class Place(), containing attri: place_id, 
                        #label, initial_tokens, tokens, continuous
                #2. collect just the values of {} self.places
                #3. for each instance in self.places.values(), take out the .tokens attri
# =============================================================================
#         print("tokens (in run_single_step())---------------------------------")
#         print(tokens) #jc/--> = exactly the 49-component list outputted to run()
#         print(len(self.places.values())) #jc/--> 49
#         print(type(self.places.values())) #jc/--> dict values
#         print(list(self.places.values())) #jc/--> list of instances of class place()
#         print(list(self.places.values())[10]) #jc/--> '24OHchol extra has 1337990.6330949916 tokens' [= the place() __str__]
#         print(type(list(self.places.values())[10])) #jc/--> 'stochastic_hfpn.Place'
#         print(list(self.places.values())[10].__dict__) #jc/--> <class 'stochastic_hfpn.Place'>
#                                                         {'place_id': 'p_24OHchol_extra', 
#                                                         'label': '24OHchol extra', 
#                                                         'initial_tokens': 4650.0, 
#                                                         'tokens': 1714948.4445276107, 
#                                                         'continuous': True}
# =============================================================================

        #Store list_of_tokens_transferred_for_each_transition_per_timestep
        #tokens_transitions = [t.list_of_produced_tokens for t in self.transitions.values()]
        list_of_tokens_transferred_for_each_transition_per_timestep_prod =[]
        list_of_tokens_transferred_for_each_transition_per_timestep_cons=[]
        #BSL: This is how you maintain the order for the panda dataframe
        for t in self.transitions.values():
            #.values() = just the values (just the instances) of self.transitions
            #--> t = instance of DiscreteTransition/ContinuousTransition()
            # print("=======creating list used in for loop") #jc/CS_CLOCK
            # print(t.transition_id)
            # print(t.list_of_consumed_tokens)
            # print("========================================")
            list_of_tokens_transferred_for_each_transition_per_timestep_prod.append(t.list_of_produced_tokens)
            list_of_tokens_transferred_for_each_transition_per_timestep_cons.append(t.list_of_consumed_tokens)
            
 
        
        # Store cumulative number of firings for each transition at specific time step
        firings = [t.firings for t in self.transitions.values()]

        # Check if time step has resulted in negative token count and raise value error 
        if any(token < 0 for token in tokens):
            # Array and list necessary to select using truth values
            place_ids = np.array(list(self.places.keys()))
            neg_token_truth_value = [token < 0 for token in tokens]
            neg_place_ids = place_ids[neg_token_truth_value]
            neg_place_tokens = np.array(tokens)[neg_token_truth_value]
            print(f"Warning: negative token count of {neg_place_tokens} in {neg_place_ids}.") #BSL:Temporarily suppress warnings
            self.counter += 1
            print(self.counter, "counter")
        for t in random_order_transitions:
            t.reset_list_of_consumed_tokens()
            t.reset_list_of_produced_tokens()
        return tokens, firings, list_of_tokens_transferred_for_each_transition_per_timestep_prod, list_of_tokens_transferred_for_each_transition_per_timestep_cons
        

    def run(self, number_time_steps, GUI_App, storage_interval=1, updategui_interval=1000):
        """Execute a run with a set amount of time-steps.

            Args:
                number_time_steps (int): number of time-steps
                storage_interval (int): number of time-steps between tokens being stored, default = 1  
                                        (e.g. storage_interval = 2: tokens stored for every 2nd time-step)
            Returns:
                single_run_tokens: 2D numpy array where first dimension is time step and second dimension is places (Row = timestep, Column = places)
                single_run_total_firings: 1D numpy array where dimension is transitions.
                df_for_rate_analytics_prod:
                df_for_rate_analytics_cons:
                FINAL_Time_step:
        """
        
        self.Dont_Save = False
        
        #jc/create 2D empty array
        single_run_tokens = np.zeros((int(number_time_steps/storage_interval)+1, len(self.places)))
        #populating the columns in the first row with the place tokens
        single_run_tokens[0] = [place.tokens for place in self.places.values()] # add initial conditions
            #2Darray[0] = calls up the 1st row of array
                
        single_run_firings = np.zeros((int(number_time_steps/storage_interval)+1, len(self.transitions)))
        single_run_firings[0] = [trans.firings for trans in self.transitions.values()]
        
        ##Make dataframe for consumption rate analytics
        list_of_column_names_consumption = []
        for trans in self.transitions.values():
            #print(global_trans_index, trans.transition_id, trans.production_coefficients) #debugging
            for input_place_id in trans.input_place_ids:
                if trans.collect_rate_analytics[0]=="yes":
                    list_of_column_names_consumption.append(trans.transition_id+" "+ input_place_id)
        
        
        
        df_for_rate_analytics_cons = pd.DataFrame(columns=list_of_column_names_consumption, index=np.arange(0,number_time_steps+1))        

        
        header_button_consump = tk.Button(GUI_App.PD_frame_in_rate_canvas, text="consumption Rates")
        header_button_consump.grid(row=0,column=0, pady=10, padx=10)

        header_button_produc = tk.Button(GUI_App.PD_frame_in_rate_canvas, text="Production Rates")
        header_button_produc.grid(row=0,column=2, pady=10, padx=10)
        
        for index, col in enumerate(df_for_rate_analytics_cons.columns):
            x = tk.Button(GUI_App.PD_frame_in_rate_canvas, text=col, state=tk.DISABLED)
            x.grid(row=index+1, column=0, padx=10, pady=10)
        
        ##Make dataframe for production rate analytics here
        list_of_column_names_production = []
        transition_count = 0 
        for global_trans_index, trans in enumerate(self.transitions.values()):
            #print(global_trans_index, trans.transition_id, trans.production_coefficients) #debugging
            for specific_trans_index, output_place_id in enumerate(trans.output_place_ids):
                transition_count = transition_count + 1
                if trans.collect_rate_analytics[1]=="yes":
                    list_of_column_names_production.append(trans.transition_id+" "+ output_place_id)
        
 
        df_for_rate_analytics_prod = pd.DataFrame(columns=list_of_column_names_production, index=np.arange(0,number_time_steps+1))
        
        for index, col in enumerate(df_for_rate_analytics_prod.columns):
            x = tk.Button(GUI_App.PD_frame_in_rate_canvas, text=col, state=tk.DISABLED)
            x.grid(row=index+1, column=2, padx=10, pady=10)      
            
            
        #GUI *****Redefine GUI Frames*****        
        
        #GUI *****Make buttons for place names*****
        list_of_places_names1 = [] 
        for key,value in self.places.items():
            list_of_places_names1.append(key)
        
        
        for index, place in enumerate(list_of_places_names1):
            tk.Button(GUI_App.frame_in_canvas, text=place, state=tk.DISABLED).grid(row=index+1, column=0, pady=10, padx=10)
            
        self.Initiate_Safe_Exit = False
        def Safe_Exit_Protocol():
            self.Initiate_Safe_Exit = True
            
            
            
            
        Exit_Button = tk.Button(GUI_App.frame_in_canvas, text="Safe Exit", command = partial(Safe_Exit_Protocol), bg="green", cursor="hand2")
        Exit_Button.grid(row=0, column=5, pady=10, padx=10)

        #Define Live Graph

        GUI_App.f = Figure(figsize=(5,5), dpi=100)
        GUI_App.a = GUI_App.f.add_subplot(111)
        GUI_App.Neuronal_Healthbar_canvas = FigureCanvasTkAgg(GUI_App.f, GUI_App.frame8)
        GUI_App.Neuronal_Healthbar_canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True) #I can also choose to grid it so its more compact for later, when I want to plot multiple plots. 
        GUI_App.toolbar = NavigationToolbar2Tk(GUI_App.Neuronal_Healthbar_canvas, GUI_App.frame8)
        GUI_App.toolbar.update() 
        self.live_plot = "False"
        def Live_Button_command():
            
            if self.live_plot == "False":
                self.live_plot = "True"
                GUI_App.Live_Button.config(text="Live!", bg="red")
            else:
                self.live_plot = "False"
                GUI_App.Live_Button.config(text="Live Graph", bg="green")
        GUI_App.Live_Button = tk.Button(GUI_App.frame8, text="Live Graph", command=Live_Button_command,bg="green", state=tk.DISABLED)
        GUI_App.Live_Button.pack(side=tk.TOP)
        
        def Plot_Live_Graph(index,t): #runs when you click on a number in run sHFPN
# =============================================================================
#             if single_run_tokens[0:int(t/storage_interval)+1,index].shape[0] != np.arange(0, (t+1)*self.time_step, self.time_step).shape[0]
#                 print("-----------------------")
#                 print("---y-------------")
#                 print(index)
#                 print(single_run_tokens.shape)
#                 print(single_run_tokens[0:int(t/storage_interval)+1,index].shape) #=y, length=t+1
#                 print("")
#                 print("---x-------------")
#                 print(t)
#                 print(self.time_step)
#                 print(np.arange(0, t*self.time_step, self.time_step).shape) #=x, length=t
#                 print(np.arange(0, t*self.time_step, self.time_step)[0])
#                 print(np.arange(0, t*self.time_step, self.time_step)[1])
#                 print(np.arange(0, t*self.time_step, self.time_step)[-2])
#                 print(np.arange(0, t*self.time_step, self.time_step)[-1])
#                 print("")
#                 print(np.arange(0, (t+1)*self.time_step, self.time_step).shape) #=x, length=t+1
#                 print(np.arange(0, (t+1)*self.time_step, self.time_step)[0])
#                 print(np.arange(0, (t+1)*self.time_step, self.time_step)[1])
#                 print(np.arange(0, (t+1)*self.time_step, self.time_step)[-2])
#                 print(np.arange(0, (t+1)*self.time_step, self.time_step)[-1])
# =============================================================================
            
            GUI_App.Live_Button.config(state="normal")
            self.plot_index = index
            the_title = str(list_of_places_names1[index])           
            GUI_App.f.clear()
            del GUI_App.a
            GUI_App.a = GUI_App.f.add_subplot(111)
            GUI_App.a.plot(np.arange(0, (t+1)*self.time_step, self.time_step), #=x
                           single_run_tokens[0:int(t/storage_interval)+1,index]) #=y
                #plot: x = np.arange(0,t*self.time_step,self.time_step)
                        #= 0...[in steps of self.time_step]...t*self.time_step
                        #= length t
                            #arange(start, stop, step_size) = give even spaced seq betw start and stop-1
                      #y = single_run_tokens[0:int(t/storage_interval)+1,index]
                        #= 1 col, length t+1, tokens from t=0 to t=now, for 1 place only
                            #index = used to refer to 1 place in single_run_tokens (= a number betw 0 and 48 inclusive)
                            #single_run_tokens = (no.timesteps+1, no.places) = eg (20001,49)
                            #0:int(t/storage_interval)+1 = take all the token vals from t=0 to t=now
                            #[0:int(t/storage_interval)+1,index] = "" for one place only
            GUI_App.a.title.set_text(the_title)
            GUI_App.a.set_xlabel('Seconds (s)')
            GUI_App.Neuronal_Healthbar_canvas.draw_idle() #BSL: note, matplotlib is not thread safe, that's why I embed it into the GUI rather than display a separate window, solving this issue was a nightmare, so I embed the live plots into the GUI.
            GUI_App.lb.itemconfig(9,{'bg':'red'}) #sets bg to red to guide user a new plot has been displayed
            
        GUI_App.f_rates = Figure(figsize=(5,5), dpi=100)
        GUI_App.a_rates = GUI_App.f_rates.add_subplot(111)
        GUI_App.Rates_canvas = FigureCanvasTkAgg(GUI_App.f_rates, GUI_App.frame10)
        GUI_App.Rates_canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True) #I can also choose to grid it so its more compact for later, when I want to plot multiple plots. 
        GUI_App.toolbar_rates = NavigationToolbar2Tk(GUI_App.Rates_canvas, GUI_App.frame10)
        GUI_App.toolbar.update() 
        self.live_plot_rates_cons = "False"  
        def Live_Button_Rates_command_Cons():
            
            if self.live_plot_rates_cons == "False":
                self.live_plot_rates_cons = "True"
                GUI_App.Live_Button_Rates.config(text="Live!", bg="red")
            else:
                self.live_plot_rates_cons = "False"
                GUI_App.Live_Button_Rates.config(text="Live Graph", bg="green")        
        GUI_App.Live_Button_Rates = tk.Button(GUI_App.frame10, text="Live Graph", command=Live_Button_Rates_command_Cons,bg="green", state=tk.DISABLED)
        GUI_App.Live_Button_Rates.pack(side=tk.TOP)        
        def Plot_Live_Graph_rates_cons(col_name,t):#arguments need to be changed
            GUI_App.Live_Button_Rates.config(state="normal")
            self.plot_col_name = col_name
            the_title = col_name           
            GUI_App.f_rates.clear()
            del GUI_App.a_rates
            GUI_App.a_rates = GUI_App.f_rates.add_subplot(111)
            numpy_array = df_for_rate_analytics_cons[col_name].to_numpy()
            GUI_App.a_rates.plot(numpy_array)
            GUI_App.a_rates.title.set_text(the_title)
            GUI_App.Rates_canvas.draw_idle() 
            GUI_App.lb.itemconfig(10,{'bg':'red'})      
            
        def Plot_Live_Graph_rates_prod(col_name,t):#arguments need to be changed
            GUI_App.Live_Button_Rates.config(state="normal")
            self.plot_col_name = col_name
            the_title = col_name           
            GUI_App.f_rates.clear()
            del GUI_App.a_rates
            GUI_App.a_rates = GUI_App.f_rates.add_subplot(111)
            numpy_array = df_for_rate_analytics_prod[col_name].to_numpy()
            GUI_App.a_rates.plot(numpy_array)
            GUI_App.a_rates.title.set_text(the_title)
            GUI_App.Rates_canvas.draw_idle() 
            GUI_App.lb.itemconfig(10,{'bg':'red'})  
            
        #GUI *****Make token buttons***** jc/on the 'run HFPN' page
        token_buttons_dict = {}
        for index,value in enumerate(single_run_tokens[0]):
            #dict keys should be the index
            index_str = str(index)
            if value > 1e6:
                readable_value = "{:e}".format(value)
            else:
                readable_value = value
            token_buttons_dict[index_str]=tk.Button(GUI_App.frame_in_canvas, text=str(readable_value), command=partial(Plot_Live_Graph,index,0), cursor="hand2")
            token_buttons_dict[index_str].grid(row=index+1, column=1, pady=10,padx=10)
        tokens_header_button = tk.Button(GUI_App.frame_in_canvas, text = "0", state=tk.DISABLED)
        tokens_header_button.grid(row=0, column=1)
        place_header_button = tk.Button(GUI_App.frame_in_canvas, text="place", state=tk.DISABLED).grid(row=0, column=0)
        GUI_App.canvas.configure(scrollregion=GUI_App.canvas.bbox("all"))
  
        
        #GUI Make Rate Token Buttons
        token_consume_rate_buttons_dict = {}
        token_produce_rate_buttons_dict = {}
        

        for index,value in enumerate(df_for_rate_analytics_cons):
            #dict keys should be the index            
            readable_value = df_for_rate_analytics_cons.iloc[0][value] #value is the column name in this case.
            index_str = str(index)
            token_consume_rate_buttons_dict[index_str]=tk.Button(GUI_App.PD_frame_in_rate_canvas, text=str(readable_value), cursor="hand2")
            token_consume_rate_buttons_dict[index_str].grid(row=index+1, column=1, pady=10,padx=10)           
        

        for index,value in enumerate(df_for_rate_analytics_prod):
            #dict keys should be the index
            readable_value = df_for_rate_analytics_prod.iloc[0][value] #value is the column name in this case.
            index_str = str(index)
            token_produce_rate_buttons_dict[index_str]=tk.Button(GUI_App.PD_frame_in_rate_canvas, text=str(readable_value), cursor="hand2")
            token_produce_rate_buttons_dict[index_str].grid(row=index+1, column=3, pady=10,padx=10)        
            
        
        # *Pause_button
        
        self.Pause_value = "off"
        def pause():
            if self.Pause_value == "off":
                Pause_button.config(text="Paused!", bg="red")
                Pause_button2.config(text="Paused!", bg="red")
                self.Pause_value = "on"
                Pickling_Button.config(state=tk.DISABLED)   
                Exit_Button.config(state=tk.DISABLED)        
            else:
                Pause_button.config(text="Pause", bg="green")
                Pause_button2.config(text="Pause", bg="green")
                self.Pause_value = "off"
                print("Unpaused!")
                Pickling_Button.config(state="normal")   
                Exit_Button.config(state="normal")
        Pause_button = tk.Button(GUI_App.frame_in_canvas, text="Pause", command = partial(pause), bg="green", cursor="hand2")
        Pause_button.grid(row=0, column=3,pady=10,padx=10)
        Pause_button2 = tk.Button(GUI_App.frame8, text="Pause", command = partial(pause), bg="green", cursor="hand2")
        Pause_button2.pack()
        

    
       
        #Premature Pickling Button
        self.Pickle_Decision = "no"
        def premature_pickle():
            self.Pickle_Decision = "yes"
        
        Pickling_Button = tk.Button(GUI_App.frame_in_canvas, text="Prematurely Pickle", command = partial(premature_pickle), cursor="hand2")
        Pickling_Button.grid(row=0, column=4,pady=10,padx=10)
        
        
#jc/ run(): firing transitions-------------------------------------------------------
    #is the actual loop: updates every time GUI updates on 'run HFPN' page
    
        for t in range(number_time_steps):
            #jc/ FIRE TRANSITIONS run_single_step() ------------------------------------
            tokens, firings,list_of_tokens_transferred_for_each_transition_per_timestep_prod, list_of_tokens_transferred_for_each_transition_per_timestep_cons  = self.run_single_step()
# =============================================================================
#             print(t) ##jc/--> incr in increments of 1, GUI updates when t reaches multiple of 1000
#             print(len(tokens)) ##jc/--> always 49
#             print("tokens (in run())--------------------------------------")
#             print(tokens) ##jc/--> list of 49 token no.s
# =============================================================================
            # print("")
            # print("======== for t in range(number_time_steps)")
            # print("--- t = " + str(t))
            
            # store current no.tokens/place ('tokens') at regular intervals (every t)            
            if t % storage_interval == 0: #jc/if current timestep is one where no. tokens should be stored into single_run_tokens (currently = EVERY t), run this
                #BSL: produce a flat list to input into pandas dataframe for later.
                flattened_list_prod = [item for sublist in list_of_tokens_transferred_for_each_transition_per_timestep_prod for item in sublist]   
                flattened_list_cons = [item for sublist in list_of_tokens_transferred_for_each_transition_per_timestep_cons for item in sublist]                  
                # print("======== create for loop") #jc/CS_CLOCK
                # print(list_of_tokens_transferred_for_each_transition_per_timestep_cons)
                
                single_run_tokens[int(t/storage_interval)+1] = tokens
                    #store tokens into the NEXT row
                        #if storage_interval=100, when t=100, want to store into 2nd row
                        #100/100 + 1 = 2
                        #200/100 + 1 = 3 = 3rd row = :)
                    #t/storage_interval = t, if storage_interval=1
                #print(single_run_tokens.shape) ##jc/--> always (total no. timesteps, no. places) = (6000001, 49)
                #print(single_run_tokens[int(t/storage_interval)+1]) ##jc/
                # print("-------- t = " + str(t)) #jc/CS_CLOCK
                # print("-------- df_for_rate_analytics_cons")
                # print(df_for_rate_analytics_cons)
                # print("")
                # print("-------- flattened_list_cons")
                # print(len(flattened_list_cons))
                # print(flattened_list_cons)
                # print("")
                # print("-------- df_for_rate_analytics_cons.loc[t+1:t+1,:]")
                # print(df_for_rate_analytics_cons.loc[t+1:t+1,:]) #whole row
                # print("")
                
                single_run_firings[int(t/storage_interval)+1] = firings
                df_for_rate_analytics_prod.loc[t+1:t+1,:] = flattened_list_prod
                df_for_rate_analytics_cons.loc[t+1:t+1,:] = flattened_list_cons
                    #df_for_rate_analytics_prod/cons = pandas dataframe
                        #cols = 1 for each arc im calculating transition for (for analytics)
                        #rows = 1 for each timestep (10k timesteps = 10k+1 rows)
                            #NB time=0 is counted as a row, with NaNs for all values
                    #flattened_list_cons/prod = list, of the rate value for each arc 
                        #in df_for_rate_analytics_prod/cons in CURRENT timestep
                    #--> sub in rate info of each arc in CURRENT timestep into the t+1th row
                        #in df_for_rate_analytics_cons dataframe
                # print("========= t % storage_interval")
                # print("-------- t = " + str(t))
                # print("single_run_tokens.shape = " + str(single_run_tokens.shape))
                
                if t % updategui_interval == 0:
                    if self.live_plot == "True":
                        try:
                            Plot_Live_Graph(self.plot_index, t)  
                        except:
                            print("NUMBER ERROR")
                            pass
                                
            # ***** GUI WIP ***** ===================
            ##Update GUI every 1k steps
            if t % updategui_interval ==0: #if t = multiple of 1000

                #Check if Paused:            
                while self.Pause_value == "on":#XXX IDK why, when GUI runs this line of code keeps running.
                    time.sleep(1)        
                    if self.Pickle_Decision=="yes":
                        Pause_button.config(state=tk.DISABLED)
                        Pause_button2.config(state=tk.DISABLED)
                        Pickling_Button.config(state=tk.DISABLED)     
                        FINAL_Time_step = t
                        break      

                tokens_header_button.config(text="Timestep: " + str(t))
                
                # print("========= t % 1000")
                # print("t = " + str(t)) 
                    #= 1000
                # print("single_run_tokens.shape = " + str(single_run_tokens.shape)) 
                    #= (101, 49) as we're creating 1 row for every storage interval
                    #= every 100th t out of 10,000 t's = 10,000/100 = 100 rows
                    #(remember we created all rows, of 0s, at beg))

                for index,value in enumerate(single_run_tokens[int(t/storage_interval)+1]):
                    #single_run_tokens = (total no. timesteps, no. places) = (10001, 49)
                    #single_run_tokens[t] = list of 49 token no.s, 1 for every place
                    #if storage_interval=100:
                        #initial condiitons --> 1st row
                        #t=0 = 1st iteration of firings have passed --> info corresponds to 2nd row in single_run_tokens
                        #t=100 --> info corresponds to 3rd row
                        #t=1000 --> 12th row 
                        #--> want single_run_tokens[11] = 1000/100 + 1
                    
                    index_str = str(index)
                    #significant_digits = 4 #for when i implement customization on the GUI to display SFs
                    if value > 1e6:
                        #value = round(value,significant_digits- int(math.floor(math.log10(abs(value)))) - 1)
                        readable_value = "{:.3e}".format(value)
                    else:
                        readable_value = value
                    token_buttons_dict[index_str].config(text=str(readable_value), command=partial(Plot_Live_Graph,index,t))
                    #Button(second_frame, text=str(value)).grid(row=index+1, column=1, pady=10,padx=10)
                    
                #Update Live Graph every 1k steps. Works but I currently disable it - update is too slow. Can probably copy and paste more to force it to update faster.
                # t_span = list(range(0,t+1,100))
                # GUI_App.a.plot(single_run_tokens[:,2])
                
                for index,col_name in enumerate(df_for_rate_analytics_cons):
                    #dict keys should be the index                    
                    value = df_for_rate_analytics_cons.iloc[t][col_name] 
                    index_str = str(index)
                    if value > 1e6:
                        readable_value = "{:e}".format(value)
                    else:
                        readable_value = value
                    token_consume_rate_buttons_dict[index_str].config(text=str(readable_value), command=partial(Plot_Live_Graph_rates_cons, col_name, t))          
                
        
                for index,col_name in enumerate(df_for_rate_analytics_prod):
                    #dict keys should be the index
                    value = df_for_rate_analytics_prod.iloc[t][col_name] 
                    index_str = str(index)
                    if value > 1e6:
                        readable_value = "{:e}".format(value)
                    else:
                        readable_value = value
                    token_produce_rate_buttons_dict[index_str].config(text=str(readable_value), command=partial(Plot_Live_Graph_rates_prod, col_name, t))
                if self.live_plot_rates_cons == "True":
                    Plot_Live_Graph_rates_cons(self.plot_col_name,t)
           
            
            #*Show final tokens in GUI when run ends --------------------------        

            if t == number_time_steps-1:
                # When t reaches number of Timestep, safe exit requirement is turned off.
                print("Debug: Turned Off Safe Exit Requirement")
                GUI_App.Safe_Exit_Required = False
                #we want to delete safe exit button
                Exit_Button.config(state=tk.DISABLED)
                # Perform Cleanup
                try:
                    GUI_App.f.clear()
                    plt.close(GUI_App.f)
                    print("Success Clear and Close F")
                except:
                    print("Unsuccessful Clear and Close F")
                    pass
                try:               
                    GUI_App.f_rates.clear()
                    plt.close(GUI_App.f_rates)
                    print("Success Clear and Close F_rates")
                except:
                    print("Unsuccessful Clear and Close F_rates")
                    pass  
                
                
                tokens_header_button.config(text="Timestep: " + str(t))
                for index,value in enumerate(single_run_tokens[int(t/storage_interval)+1]):
                    index_str = str(index)
                    #significant_digits = 4 #for when i implement customization on the GUI to display SFs
                    if value > 1e6:
                        #value = round(value,significant_digits- int(math.floor(math.log10(abs(value)))) - 1)
                        readable_value = "{:.3e}".format(value)
                    else:
                        readable_value = value
                    token_buttons_dict[index_str].config(text=str(readable_value))            
                
                self.Pause_value =="off"
                Pause_button.config(state=tk.DISABLED)
                Pause_button2.config(state=tk.DISABLED)
                Pickling_Button.config(state=tk.DISABLED)
                GUI_App.Live_Button.config(state=tk.DISABLED)
                FINAL_Time_step = t
                    

            #Premature Pickle
            if self.Pickle_Decision =="yes":
                Pause_button.config(state=tk.DISABLED)
                Pause_button2.config(state=tk.DISABLED)
                Pickling_Button.config(state=tk.DISABLED)
                GUI_App.Live_Button.config(state=tk.DISABLED)
                FINAL_Time_step = t
                Exit_Button.destroy()
                
                try:
                    GUI_App.f.clear()
                    plt.close(GUI_App.f)
                    print("Success Clear and Close F")
                except:
                    print("Unsuccessful Clear and Close F")
                    pass
                try:               
                    GUI_App.f_rates.clear()
                    plt.close(GUI_App.f_rates)
                    print("Success Clear and Close F_rates")
                except:
                    print("Unsuccessful Clear and Close F_rates")
                    pass 
            
                
                break
            if self.Initiate_Safe_Exit == True:
                try:
                    GUI_App.f.clear()
                    plt.close(GUI_App.f)
                    print("Success Clear and Close F")
                except:
                    print("Unsuccessful Clear and Close F")
                    pass
                try:               
                    GUI_App.f_rates.clear()
                    plt.close(GUI_App.f_rates)
                    print("Success Clear and Close F_rates")
                except:
                    print("Unsuccessful Clear and Close F_rates")
                    pass  
                FINAL_Time_step = t
                Exit_Button.destroy()
                Pause_button.config(state=tk.DISABLED)
                Pickling_Button.config(state=tk.DISABLED)
                
                
                break
        if self.Initiate_Safe_Exit == True:
            print("Initiated Safe Exit so no Saving")
            self.Dont_Save = True
        
                
        #Disable graphing ability after running so user doesnt create embedded graphs which cause errors later.
        for index,value in enumerate(single_run_tokens[0]):
            #dict keys should be the index
            index_str = str(index)
            if value > 1e6:
                readable_value = "{:e}".format(value)
            else:
                readable_value = value
            token_buttons_dict[index_str].config(state=tk.DISABLED) 

        
        for index,col_name in enumerate(df_for_rate_analytics_cons):
            #dict keys should be the index                    
            value = df_for_rate_analytics_cons.iloc[t][col_name] 
            index_str = str(index)
            if value > 1e6:
                readable_value = "{:e}".format(value)
            else:
                readable_value = value
            token_consume_rate_buttons_dict[index_str].config(state=tk.DISABLED)      
        

        for index,col_name in enumerate(df_for_rate_analytics_prod):
            #dict keys should be the index
            value = df_for_rate_analytics_prod.iloc[t][col_name] 
            index_str = str(index)
            if value > 1e6:
                readable_value = "{:e}".format(value)
            else:
                readable_value = value
            token_produce_rate_buttons_dict[index_str].config(state=tk.DISABLED)

        # Determine how many times transition fired in single run
        single_run_total_firings = single_run_firings[-1,:]
       
        GUI_App.lb.itemconfig(7, fg="lime")
        return single_run_tokens, single_run_total_firings, df_for_rate_analytics_prod,  df_for_rate_analytics_cons, FINAL_Time_step



    def run_many_times(self, number_runs, number_time_steps, GUI_App, storage_interval=1, updategui_interval=1000):
        """Runs multiple iterations of the HFPN.
        
            Args: 
                number_runs (int): total number of runs
                number_time_steps (int): number of time steps for each run
                storage_interval (int): number of time-steps between tokens being stored, default = 1  
                                        (e.g. storage_interval = 2: tokens stored for every 2nd time-step)
        """
                
        if storage_interval == -1: #!!!
            storage_interval = max(1, int(number_time_steps/1000))
        
        if number_time_steps % storage_interval != 0:
            raise ValueError('Number of time-steps should be a multiple of 1000')

        # Store time (in seconds) at each time step 
        self.time_array = np.arange(0, self.time_step*(number_time_steps+1), self.time_step*storage_interval)
        #time_step*number_time_steps gives you the total time.
        #print(self.time_array) brandoggy
        
        # First dimension = run number, second dimension = time step, third dimension = place 
        self.token_storage = np.zeros((number_runs, int(number_time_steps/storage_interval)+1, len(self.places)))
            #jc/ array of 0s, first dimension = 1, since im only doing 1 run
        #print(self.token_storage) ##jc/test

        # First dimension = run number, second dimension = transition
        self.firings_storage = np.zeros((number_runs, len(self.transitions)))
     
        
# jc/ run_many_times: run() ------------------------------------------------------------------       
#storing it for visualisation
        for i in range(number_runs):
            #jc/ will loop only once thru for me, if only doing 1 run
            self.token_storage[i], self.firings_storage[i], self.df_for_rate_analytics_prod, self.df_for_rate_analytics_cons, self.FINAL_time_step = self.run(number_time_steps, GUI_App, storage_interval, updategui_interval)
                #jc/ method self.run() returns 5 outputs. 1st two are stored inside the
                    #2 arrays of zeroes (self.token_storage, etc)
                    #self.firings_storage = array of 0s and 1s, for each place (col) for each timestep (row)
                        #1 = that place fires in this timestep
            #print(self.firings_storage) ##jc/test
            self.reset_network()
#-------------------------------------------------------------------------------------------

        # Store mean number of firings for each transition across all runs
        self.mean_firings = np.mean(self.firings_storage, axis = 0)
        
        #store the number of timesteps used 
        self.the_number_of_timesteps = number_time_steps
        
        

        
    def reset_network(self):
        """Resets the token values for each place and number of firings for each transition."""
        # Loop through places and reset the network 
        for place in self.places.values():
            place.reset()
        # Loop through transitions and reset number of firings
        for t in self.transitions.values():
            t.reset() 

#    def GUI_plot(place_id, token_storage):
        
            
   
# if __name__ == '__main__':

#     # Initialize an empty Petri net
#     pn = HFPN(time_step=0.001, printout=True)

#     # Add places for each chemical species
#     pn.add_place(initial_tokens=100, place_id="p_H2", label="Hydrogen", continuous=True)
#     pn.add_place(100, place_id="p_O2", label="Oxygen", continuous=True)
#     pn.add_place(0, place_id="p_H2O", label="Water", continuous=True)
#     pn.add_place(0, place_id="p_I", label="Inhibitor", continuous=True)

#     rate_constant = 0.001

#     ### Add the same transition in three different ways:
#     # 1. Directly using speed_functions
#     # pn.add_transition(  transition_id = 't_a',
#     #                     label = 'Example transition a',
#     #                     input_place_ids = ['p_H2', 'p_O2', 'p_I'],
#     #                     firing_condition = lambda a : (a['p_H2'] >= 0 or a['p_O2'] >= 0) and a['p_I'] <= 0.01,
#     #                     consumption_speed_functions = [lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 2, 
#     #                                                    lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 1,
#     #                                                    lambda a : 0],
#     #                     output_place_ids = ['p_H2O'],
#     #                     production_speed_functions = [lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 2],
#     #                     delay = 0.00
#     # )

#     # 2. Using one shared reaction_speed_function for each species 
#     pn.add_transition_with_speed_function(
#                         transition_id = 't_b',
#                         label = 'Example transition b',
#                         input_place_ids = ['p_H2', 'p_O2', 'p_I'],
#                         firing_condition = lambda a : a['p_H2'] >= 0 or a['p_O2'] >= 0 and a['p_I'] <= 0.01,
#                         reaction_speed_function = lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1,
#                         consumption_coefficients = [20, 10, 0], 
#                         output_place_ids = ['p_H2O'],
#                         production_coefficients = [2],
#                         stochastic_parameters=[1])

#     # # 3. Using mass-action as the shared reaction_speed_function
#     # pn.add_transition_with_mass_action(  transition_id = 't_c',
#     #                     label = 'Example transition c',
#     #                     rate_constant = rate_constant,
#     #                     input_place_ids = ['p_H2', 'p_O2', 'p_I'],
#     #                     firing_condition = lambda a : a['p_H2'] >= 0 and a['p_O2'] >= 0 and a['p_I'] <= 0.01,
#     #                     consumption_coefficients = [2, 1, 0],
#     #                     output_place_ids = ['p_H2O'],
#     #                     production_coefficients = [2]
#     # )

#     # # Adding transition using Michaelis Menten
#     # pn.add_transition_with_michaelis_menten(transition_id = 't_michaelis_menten',
#     #                                     label = 'Michaelis Menten test',
#     #                                     Km = 1,
#     #                                     vmax = 1,
#     #                                     input_place_ids = ['p_H2', 'p_O2'],
#     #                                     substrate_id = 'p_H2',
#     #                                     consumption_coefficients = [1, 0],
#     #                                     output_place_ids = ['p_H2O'],
#     #                                     production_coefficients = [1],
#     #                                     vmax_scaling_function = lambda a : 1)


#     pn.run_many_times(1,5, -1)
#    # print(pn.transitions['t_b'].label)
#     #print('shape of token_storage:', pn.token_storage.shape)
#     #print('time array:', pn.time_array)
#     #print(pn.token_storage)
#     print(pn.list_of_list_of_tokens_transferred_for_each_transition_per_timestep)
