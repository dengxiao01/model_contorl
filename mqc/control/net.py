# -*- coding: utf-8 -*-

import re 
from mqc.utils import *
from mqc.defaults import * 

class Nets():
    """
    Get net substance output information.

    """
    def __init__(self):
        """
        Define net substance output.

        """
        self.nets = []

    def set_no_c_source_supply(model_info, model):
        """
        Set no C source supply, change the reaction boundary outside the model exchange reaction to (-1000, 1000) within this range.

        """
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met = rxn["all_mets"][0]
                if any(r_met.startswith(source_name) for source_name in AUTOTROPHIC_SOURCE):
                    model.reactions.get_by_id("rxn['id']").bounds = (0,1000)
                    rxn["bounds"] = [0,1000]
                    rxn["autotrophic"] = "true"
                if re.search("C[A-Z]",r_met.formula) or re.search("C[\d]",r_met.formula) or r_met.formula.find('R') != -1 or r_met.formula.find('X') != -1: # C后边是大写字母或者数字
                    if rxn["id"] not in CO2 or rxn["id"] not in HCO3:
                        model.reactions.get_by_id("rxn['id']").bounds = (0,1000)
                        rxn["bounds"] = [0,1000]
                        rxn["carbon_source"] = "true"
            else:
                if rxn["bounds"][0] < 0 and rxn["bounds"][1] > 0:
                    model.reactions.get_by_id("rxn['id']").bounds = (-1000,1000)
                    rxn["bounds"] = [-1000,1000]
                if rxn["bounds"][0] > 0 and rxn["bounds"][1] > 0:
                    model.reactions.get_by_id("rxn['id']").bounds = (0,1000)
                    rxn["bounds"] = [0,1000]
                if rxn["bounds"][0] < 0 and rxn["bounds"][1] < 0:
                    model.reactions.get_by_id("rxn['id']").bounds = (-1000,0)
                    rxn["bounds"] = [-1000,0]


    def get_atpm_rxn(self, model_info, model):
        """
        get atpm reactions
        """
        if set(model_info["all_rxn_obj"]) & set(ATPM) == 0:
            atpm_id = add_atpm_rxn(model)
        else:
            atpm_id = [i for i in ATPM if i in model_info["all_rxn_obj"]][0]
        return atpm_id


    def find_net_generation(model_info, model):
        r_met, net_generation, met_num, scoring = [], [], 0, 0
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"][0]
        for met in model['metabolites']:
            met_num += 1
            if met['name'] not in FREE_METS:
                if re.search("C[A-Z]",met['formula']) or re.search("C[\d]",met['formula']): # C后边是大写字母或者数字
                    objectiveId = add_demand(model_info, model, met['id'])
                    with model:
                        if 'DM_' + str(met.id) in model.reactions:
                            model.objective = 'DM_' + str(met.id)
                            solution = model.optimize()
                            if solution.fluxes['DM_' + str(met.id)] > 1e-5 and '_'.join(str(met).split('_')[:-1]) not in r_met:
                                net_generation.append(met.id)
                                scoring += 1
                        else:
                            demand_rxn = model.add_boundary(met, type='demand') # 添加demand反应
                            model.objective = demand_rxn.id
                            solution = model.optimize()
                            if solution.fluxes[demand_rxn.id] > 1e-5 and '_'.join(str(met).split('_')[:-1]) not in r_met:#操作同上
                                net_generation.append(met.id)
                                scoring += 1
        scoring = (met_num - scoring) / met_num



    def net_control(self, model_info, model, check_model, model_control_info):
        """
        nets correction process
        """
        if not is_bio_exp_correct(model_info):
            model_control_info["biomass_info"] = "无法找到biomass方程,请检查后再试!!!"
            return 0
        self.set_no_c_source_supply(model_info, model)
        atpm_id = self.get_atpm_rxn()
        set_model_objective(model_info, model, atpm_id)