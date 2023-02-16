# -*- coding: utf-8 -*-

from mqc.control.model_control import ModelPreprocess
from mqc.control.preprocessing_control import Preprocess
from mqc.utils import *

class Biomass():


    def __init__(self):
        """
        Get model information and check_model object.

        """
        # modelpreprocess = ModelPreprocess()
        # modelpreprocess.get_model_info()
        # self.modelInfo = modelpreprocess.model_info
        # self.check_model = modelpreprocess.check_model
        self.biomass_info = []

    


    def fix_biomass(self, model_info, model, check_model, model_control_info):
        """
        Biomass correction process
        """
        if not is_bio_exp_correct(model_info):
            model_control_info["biomass_info"] = "无法找到biomass方程,请检查后再试!!!"
            return 0
        
