# -*- coding: utf-8 -*-

import os, sys
import json

from utils import *
# from control import preprocessing_control, model_control, net_control
sys.path.append('/home/dengxiao/workspace/model_contorl')

from mqc.control.preprocessing_control import Preprocess

from control.model_control import ModelPreprocess
from control.nadh_control import Nadhs
from control.atp_control import Atps
from control.net_control import Nets
from control.yield_control import Yields
from control.biomass_control import Biomasses



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
        self.model_control_info["nadhs"] = {}
        self.model_control_info["atps"] = {}
        self.model_control_info["nets"] = {}    
        self.model_control_info["yields"] = {}
        self.model_control_info["biomasses"] = {}


    def write_result(self, model_control_info):
        """"""
        result = json.dumps(model_control_info, ensure_ascii=False, allow_nan = True, indent=1)
        with open("mqc/result.json", "w", newline='\n',) as f:
            f.write(result)

    def write_result2(self, model_control_info, model):
        """"""
        result = json.dumps(model_control_info, ensure_ascii=False, allow_nan = True, indent=1)
        with open(f"mqc/result/{model}.json", "w", newline='\n',) as f:
            f.write(result)


    def model_control(self):
        """
        The overall process of model quality control.

        """
        file_path = get_model_file()
        controler = Preprocess(file_path)
        model_pre_process = ModelPreprocess()
        model_pre_process.get_model_info(controler)

        nadhs = Nadhs()
        if nadhs.nadh_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) == 0:
            return 0
        self.write_result(self.model_control_info)
        atps = Atps()
        atps.atp_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) == 0
        self.write_result(self.model_control_info)
        nets = Nets()
        nets.net_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)
        self.write_result(self.model_control_info)
        yields = Yields()
        yields.yield_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) 
        self.write_result(self.model_control_info)
        biomasses = Biomasses()
        biomasses.biomass_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) 
  

        # with open("mqc/result.json", "w") as f:
        #     json.dump(self.model_control_info, f, ensure_ascii=False)
        
        self.write_result(self.model_control_info)
        
        

    def model_control2(self):
        """
        The overall process of model quality control.

        """
        model_ok = ['RECON1.xml','Recon3D.xml','e_coli_core.xml','iAB_RBC_283.xml','iAF1260.xml','STM_v1_0.xml']
        model_ok2 = os.listdir("mqc/result")
        files = os.listdir("mqc/test_data/bigg_data")
        try:
            for file in files:
                if file.split('.')[0]+'.json' not in model_ok2:
                    file_path = f"mqc/test_data/bigg_data/{file}"
                    controler = Preprocess(file_path)
                    model_pre_process = ModelPreprocess()
                    model_pre_process.get_model_info(controler)

                    nadhs = Nadhs()
                    if nadhs.nadh_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) == 0:
                        break 
                    self.write_result2(self.model_control_info)
                    atps = Atps()
                    atps.atp_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)
                    self.write_result2(self.model_control_info, controler.model)
                    nets = Nets()
                    nets.net_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)
                    self.write_result2(self.model_control_info, controler.model)
                    yields = Yields()
                    yields.yield_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) 
                    self.write_result2(self.model_control_info, controler.model)
                    biomasses = Biomasses()
                    biomasses.biomass_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info) 
            

                    # with open("mqc/result.json", "w") as f:
                    #     json.dump(self.model_control_info, f, ensure_ascii=False)
                    
                    # self.write_result(self.model_control_info)
                    self.write_result2(self.model_control_info, controler.model)
        except Exception as e:
            print(controler.model.id, repr(e))


modelControl = ModelControl()
modelControl.model_control()