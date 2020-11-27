import math
import simulation_settings as ss
import constants
import numpy as np

class Driver:
    shapes = []
    params = []
    num_layers = 0
    num_nodes = 0

    def __init__(self,shapes):
        self.params = []
        self.shapes = shapes
        for s in self.shapes:
            temp = np.zeros(s)
            self.params.append(temp)

        self.num_layers = shapes.shape[0]
    
    def set_param(self,new_params):
        if(params.shape == new_params.shape):
            params = new_params
    
    def output(self,input_array):
        input_layer = input_array
        temp = None
        output = 0
        for i in range(self.num_layers):
            # print(temp)
            # print(i)
            if(i==0):
                temp = np.matmul(input_layer,self.params[0])
            elif(i==self.num_layers-1):
                output = np.matmul(temp,self.params[i])
                print(output.shape)
            else:
                temp = np.matmul(temp,self.params[i])
            # print(temp)
        return output[0]