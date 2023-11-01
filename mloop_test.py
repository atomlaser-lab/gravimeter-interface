# -*- coding: utf-8 -*-
import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv

#Other imports
import numpy as np
import scipy.io
import time
import os

#Declare your custom class that inherits from the Interface class
class CustomInterface(mli.Interface):
    
    #Initialization of the interface, including this method is optional
    def __init__(self):
        #You must include the super command to call the parent class, Interface, constructor 
        super(CustomInterface,self).__init__()
        
        #Attributes of the interface can be added here
        #If you want to precalculate any variables etc. this is the place to do it
        #In this example we will just define the location of the minimum
        # self.minimum_params = np.array([0,0.1,-0.1])
        
    #You must include the get_next_cost_dict method in your class
    #this method is called whenever M-LOOP wants to run an experiment
    def get_next_cost_dict(self,params_dict):
        
        #Get parameters from the provided dictionary
        # params = params_dict['params']
        
        #Here you can include the code to run your experiment given a particular set of parameters
        #In this example we will just evaluate a sum of sinc functions
        scipy.io.savemat('exp_input.mat',params_dict)
        while not os.path.isfile('exp_output.mat'):
            time.sleep(1)
            
        mat = scipy.io.loadmat('exp_output.mat',mat_dtype=True,squeeze_me=True)
        os.remove('exp_output.mat')
        
        #The cost, uncertainty and bad boolean must all be returned as a dictionary
        #You can include other variables you want to record as well if you want
        cost_dict = {'cost':mat['cost'], 'uncer':mat['uncer'], 'bad':mat['bad']}
        return cost_dict
    
def main():
    #M-LOOP can be run with three commands
    
    #First create your interface
    interface = CustomInterface()
    #Next create the controller. Provide it with your interface and any options you want to set
    controller = mlc.create_controller(interface, 
                                       max_num_runs = 300,
                                       max_num_runs_without_better_params = 50,
                                       num_params = 12, 
                                       min_boundary = [-1,  -1, -1, -1, -1, -1, -1, -1, 0,  0,  0,  0],
                                       max_boundary = [1,   1,  1,  1,   1,  1,  1,  1,  0.1,    0.1,    0.1,    0.1],
                                       first_params = [0,   0,  0,  0,   0,  0,  0,  0,  20e-3,   20e-3,   20e-3,   20e-3],
                                       trust_region = 0.4,
                                       no_delay = False,
                                       controller_archive_file_type = 'mat',
                                       learner_archive_file_type = 'mat')
    #To run M-LOOP and find the optimal parameters just use the controller method optimize
    controller.optimize()
    
    #The results of the optimization will be saved to files and can also be accessed as attributes of the controller.
    print('Best parameters found:')
    print(controller.best_params)
    
    #You can also run the default sets of visualizations for the controller with one command
    mlv.show_all_default_visualizations(controller)
    

#Ensures main is run when this code is run as a script
if __name__ == '__main__':
    main()