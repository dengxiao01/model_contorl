# -*- coding: utf-8 -*-

import os, sys
# from mqc.control import preprocessing_control 
# sys.path.append(os.getcwd)
# from os.path import abspath, join, dirname
# sys.path.insert(0, join(abspath(dirname(__file__)), 'mqc'))
# print(os.path.abspath(__file__))
import json

from utils import *
# from control import preprocessing_control, model_control, net_control
sys.path.append('/home/dengxiao/workspace/model_contorl')
# from mqc import 
from mqc.control.preprocessing_control import Preprocess

from control.model_control import ModelPreprocess
from control.net_control import Nets


__version__ = "0.1.0"
__author__ = "dengxiao"

# PROJECT_ABSOLUTE_PATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# print('PROJECT_ABSOLUTE_PATH:%s\n' % PROJECT_ABSOLUTE_PATH)
# print(os.path.dirname(os.path.abspath(__file__)))
# sys.path.append(PROJECT_ABSOLUTE_PATH)
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

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

        with open("mqc/result.json", "w") as f:
            json.dump(self.model_control_info, f, ensure_ascii=False)

        # biomassprocess = Biomass()
        # biomassprocess.fix_biomass(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)

modelControl = ModelControl()
modelControl.model_control()