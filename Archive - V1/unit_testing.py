import unittest
import numpy as np

from utils.petri_nets import PetriNet, PetriNetModel
from utils.hfpn import HFPN

class TestHFPN(unittest.TestCase):

    def test_add_transition_methods(self):

        # Initialize an empty Petri net
        pn1 = HFPN()
        pn2 = HFPN()
        pn3 = HFPN()

        # Add places for each chemical species
        pn1.add_place(initial_tokens=20, place_id="p_H2", label="Hydrogen", continuous=True)
        pn1.add_place(20, place_id="p_O2", label="Oxygen", continuous=True)
        pn1.add_place(0, place_id="p_H2O", label="Water", continuous=True)
        pn1.add_place(0, place_id="p_I", label="Inhibitor", continuous=True)

        # Add places for each chemical species
        pn2.add_place(initial_tokens=20, place_id="p_H2", label="Hydrogen", continuous=True)
        pn2.add_place(20, place_id="p_O2", label="Oxygen", continuous=True)
        pn2.add_place(0, place_id="p_H2O", label="Water", continuous=True)
        pn2.add_place(0, place_id="p_I", label="Inhibitor", continuous=True)

        # Add places for each chemical species
        pn3.add_place(initial_tokens=20, place_id="p_H2", label="Hydrogen", continuous=True)
        pn3.add_place(20, place_id="p_O2", label="Oxygen", continuous=True)
        pn3.add_place(0, place_id="p_H2O", label="Water", continuous=True)
        pn3.add_place(0, place_id="p_I", label="Inhibitor", continuous=True)

        rate_constant = 0.3

        ### Add the same transition in three different ways:

        # 1. Directly using speed_functions
        pn1.add_transition(  transition_id = 't_a',
                            label = 'Example transition a',
                            input_place_ids = ['p_H2', 'p_O2', 'p_I'],
                            firing_condition = lambda a : a['p_H2'] >= 0 or a['p_O2'] >= 0 and a['p_I'] <= 0.01,
                            consumption_speed_functions = [lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 2, 
                                                           lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 1,
                                                           lambda a : 0],
                            output_place_ids = ['p_H2O'],
                            production_speed_functions =  [lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1 * 2]
        )

        # 2. Using one shared reaction_speed_function for each species 
        pn2.add_transition_with_speed_function(
                            transition_id = 't_b',
                            label = 'Example transition b',
                            input_place_ids = ['p_H2', 'p_O2', 'p_I'],
                            firing_condition = lambda a : a['p_H2'] >= 0 or a['p_O2'] >= 0 and a['p_I'] <= 0.01,
                            reaction_speed_function = lambda a : rate_constant * a['p_H2']**2 * a['p_O2']**1,
                            consumption_coefficients = [2, 1, 0], 
                            output_place_ids = ['p_H2O'],
                            production_coefficients = [2])

        # 3. Using mass-action as the shared reaction_speed_function
        pn3.add_transition_with_mass_action(  transition_id = 't_c',
                            label = 'Example transition c',
                            rate_constant = rate_constant,
                            input_place_ids = ['p_H2', 'p_O2', 'p_I'],
                            firing_condition = lambda a : a['p_H2'] >= 0 and a['p_O2'] >= 0 and a['p_I'] <= 0.01,
                            consumption_coefficients = [2, 1, 0],
                            output_place_ids = ['p_H2O'],
                            production_coefficients = [2]
        )

        pn1.run_many_times(number_runs=1, number_time_steps=10)
        pn2.run_many_times(number_runs=1, number_time_steps=10)
        pn3.run_many_times(number_runs=1, number_time_steps=10)

        self.assertTrue(np.array_equal(pn1.token_storage, pn2.token_storage))
        self.assertTrue(np.array_equal(pn1.token_storage, pn3.token_storage))


class TestPetriNets(unittest.TestCase):

    def test_single_simple_net(self):
        # Initialize an empty Petri net
        pn = PetriNetModel()

        # Add places for each chemical species
        pn.add_place(20, place_id="p_H2", label="hydrogen")
        pn.add_place(20, place_id="p_O2", label="oxygen")
        pn.add_place(0, place_id="p_H2O", label="water")

        # Add transition corresponding to chemical reaction
        pn.add_transition(  transition_id       = "t_water",
                            label 	            = "Transition 1",
                            input_place_ids	    = ['p_H2', 'p_O2'],
                            input_arc_weights   = [2, 1], 
                            output_place_ids    = ['p_H2O'],
                            output_arc_weights  = [2])

        # Run the network X times
        [pn.run_one_step() for _ in range(12)]

        # Check firing sequence is as expected
        self.assertEqual(pn.successful_firings, [True for _ in range(10)] + [False for _ in range(2)])

        # Check final token counts are as expected
        self.assertEqual(pn.places['p_H2'].tokens, 0)
        self.assertEqual(pn.places['p_O2'].tokens, 10)
        self.assertEqual(pn.places['p_H2O'].tokens, 20)


    # Initializes a collection of petri net copies using the main PetriNet method. 
    def test_collection_of_simple_nets(self):
        pn = PetriNet(10) # run 10 copies of net and see outcome

        # Add places for each chemical species
        pn.add_place(20, place_id="p_H2", label="hydrogen")
        pn.add_place(20, place_id="p_O2", label="oxygen")
        pn.add_place(0, place_id="p_H2O", label="water")

        # Add transition corresponding to chemical reaction
        pn.add_transition(  transition_id       = "t_water",
                            label 	            = "Transition 1",
                            input_place_ids	    = ['p_H2', 'p_O2'],
                            input_arc_weights   = [2, 1], 
                            output_place_ids    = ['p_H2O'],
                            output_arc_weights  = [2])

        # Run the network X times
        pn.run(12)

        # Check final token counts are as expected
        self.assertEqual(pn.timeseries_mean_for_place('p_H2')[-1], 0)
        self.assertEqual(pn.timeseries_mean_for_place('p_O2')[-1], 10)
        self.assertEqual(pn.timeseries_mean_for_place('p_H2O')[-1], 20)
        self.assertEqual([s for s in pn.timeseries_std[-1]], [0, 0, 0]) # std = 0 for deterministic net


if __name__ == '__main__':
    
    unittest.main()

