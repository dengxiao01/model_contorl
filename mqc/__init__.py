# -*- coding: utf-8 -*-

import sys 



__version__ = "0.1.0"
__author__ = "dengxiao"

from control import *

class ModelControl():
    """
    Obtain the total output information of the model quality control.

    """
    def __init__(self):
        """
        define output dictionary.

        """
        self.model_control_info = {}

    def model_control(self):
        """
        The overall process of model quality control.

        """
        file_path = get_model_file()
        controler = Preprocess(file_path)
        model_pre_process = ModelPreprocess()
        model_pre_process.get_model_info(controler)
        nets = Nets()
        nets.net_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)




        biomassprocess = Biomass()
        biomassprocess.fix_biomass(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)