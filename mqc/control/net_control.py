# -*- coding: utf-8 -*-

import re 
from mqc.utils import *
from mqc.defaults import * 
from mqc.control.rules_control import Rules

class Nets():
    """
    Get net substance output information.

    """
    def __init__(self):
        """
        Define net substance output.

        """
        self.nets = []

    def set_no_c_source_supply(self, model_info, model):
        """
        Set no C source supply, change the reaction boundary outside the model exchange reaction to (-1000, 1000) within this range.

        """
        for rxn in model_info["reactions"]:
            rxns = model.reactions.get_by_id(rxn['id'])
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met_name = rxn["all_mets"][0] 
                r_met = list(rxns.metabolites.keys())[0]
                if any(r_met_name.startswith(source_name) for source_name in AUTOTROPHIC_SOURCE):
                    rxns.bounds = (0,1000)
                    rxn["bounds"] = [0,1000]
                    rxn["autotrophic"] = "true"
                if re.search("C[A-Z]",r_met.formula) or re.search("C[\d]",r_met.formula) or r_met.formula.find('R') != -1 or r_met.formula.find('X') != -1: # C后边是大写字母或者数字
                    if r_met_name not in CO2_NAME or r_met_name not in HCO3_NAME:
                        rxns.bounds = (0,1000)
                        rxn["bounds"] = [0,1000]
                        rxn["carbon_source"] = "true"
            else:
                if rxn["bounds"][0] < 0 and rxn["bounds"][1] > 0:
                    rxns.bounds = (-1000,1000)
                    rxn["bounds"] = [-1000,1000]
                if rxn["bounds"][0] > 0 and rxn["bounds"][1] > 0:
                    rxns.bounds = (0,1000)
                    rxn["bounds"] = [0,1000]
                if rxn["bounds"][0] < 0 and rxn["bounds"][1] < 0:
                    rxns.bounds = (-1000,0)
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


    def find_net_generation(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        r_met, net_generation, met_num, scoring = [], [], 0, 0
        for rxn in model_info["reactions"]:  # Find the metabolites of the reactions with the lowest bound less than 0 in the exchange reaction
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"]
        for met in model_info['metabolites']:
            met_num += 1
            if met['name'] not in FREE_METS:
                if re.search("C[A-Z]",str(met['formula'])) or re.search("C[\d]",str(met['formula'])): # Carbon followed by capital letters or numbers         
                    with model:
                        objectiveId = add_demand(model_info, model, met['id'])
                        model.objective = objectiveId
                        if model.optimize().objective_value > 1e-5 and met['name'] not in r_met:
                            net_generation.append(met['id'])
                            scoring += 1
        scoring = (met_num - scoring) / met_num
        net_proportion = f'{(1 - scoring) : .2%}'
        model_control_info["nets"]["净物质分数:"] = f"{scoring :.2%}"
        model_control_info["nets"]["代谢物总数:"] = met_num
        model_control_info["nets"]["净物质总数:"] = len(net_generation)
        model_control_info["nets"]["净物质占比:"] = net_proportion
        model_control_info["nets"]["净物质:"] = net_generation


    def net_penalty_points(self, model_info, model, rules):
        """
        give reaction penalty
        """
        rules.get_all_rules(model_info, model)
        check_rxn_balance(model_info, model)
        check_C_balance(model_info, model)
        for rxn in model_info['reactions']:
            grate = 0
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns'] and not rxn['id'].startswith('ADD_'):
                if rxn['balance'] == 'false':  # Response Imbalance Response Penalty
                    grate += 3
                if rxn['c_balance'] == 'false':  # carbon imbalance reaction penalty
                    grate += 20
                if rxn["rules"]:  # 3 points for failure to comply with reaction direction rules
                    grate += 3
                rxn["net_penalty_points"] = grate
            else:
                rxn["net_penalty_points"] = 0


    def close_max_value_net(self, model_info, model_control_info, model, object_rxn_id, netId):
        """
        Turn off the maximum value in each result
        """
        grate_delete_rxn = []
        count = 0
        if model.optimize().fluxes[object_rxn_id] <= 1e-5:
            print(model,'.......',object_rxn_id,'..........yes')
        while model.optimize().fluxes[object_rxn_id] > 1e-5:  
            count += 1
            print(count,' -- ',object_rxn_id)
            if count > 100:
                print('error: 净物质无限循环！')
                model_control_info["net"] = "净物质无限循环！"
                exit()
            write_flux_file(model_info, model)
            fba_rxn_dic={}
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            for ids,v in need_fluxes.items():      # i : 'ATPM' 'ADD_ATPM'
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        rxns = model.reactions.get_by_id(ids)
                        reactants_mets=[m.id for m in rxns.reactants]
                        products_mets=[m.id for m in rxns.products]
                        if rxn['balance'] == "false":  # 生成目标物质的反应不能关闭
                            if (v > 0 and netId in products_mets) or (v < 0 and netId in reactants_mets):
                                continue 
                        if ids == object_rxn_id: # 当前结果文件里的DM目标反应不能参与下面的关闭
                            continue
                        fba_rxn_dic[ids] = rxn['net_penalty_points']  # 'ATPM':3  'ADD_ATPM':0

            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()):#关掉每次结果中的所有最大值
                    if ids != 'ADD_ATPM' and ids != 'ATPM' and fba_rxn_dic[ids]!= 0:  # 目标反应不能关
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(model.reactions.get_by_id(ids),fba_rxn_dic[ids])
            print('close :',model.optimize().fluxes[object_rxn_id])
        return grate_delete_rxn


    def add_back_reactions_net(self, model_info, model, check_model, grate_delete_rxn):
        """
        Add the checked reaction back to see if NADH can still be generated infinitely
        """
        print('grate_delete_rxn: ',grate_delete_rxn)
        infinite_rxn = []
        for rxnId in grate_delete_rxn:
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(rxnId).bounds
            if model.optimize().objective_value >= 1e-6:
                pfba_solution=pfba(model)
                rxns = model.reactions.get_by_id(rxnId)  
                for rxn in model_info['reactions']:
                    rxn['net_modify'] = 'false'
                    if rxnId == rxn['id']:
                        if rxn['balance'] == 'false':
                            rxns.bounds = (0,0)
                            rxn['bounds'] = [0,0]
                            rxn['net_modify'] = "true"
                            infinite_rxn.append(rxnId)
                        else:
                            if pfba_solution.fluxes[rxnId] > 0:
                                rxns.bounds = (-1000,0)
                                rxn['bounds'] = [-1000,0]
                                rxn['net_modify'] = "true"
                                infinite_rxn.append(rxnId)
                            if pfba_solution.fluxes[rxnId] < 0:
                                rxns.bounds = (0,1000)
                                rxn['bounds'] = [0,1000]
                                rxn['net_modify'] = "true"
                                infinite_rxn.append(rxnId)
        return infinite_rxn


    def get_net_infinite_info(self, model_info, model, model_control_info):
        """
        Get response corrections
        """
        infinite_info, bounds_modification = [], ''
        for rxn in model_info['reactions']:
            if rxn['net_modify'] == 'true' and rxn['rules']:
                rxns = model.reactions.get_by_id(rxn['id'])
                if rxn['bounds'] == [0,0] : bounds_modification = "close"
                if rxn['bounds'] == [-1000,0] : bounds_modification = "<--"
                if rxn['bounds'] == [0,1000] : bounds_modification = "-->"
                infiniteInfo = {"反应ID" : rxn["id"],
                          "反应表达式" : rxn['rxn_exp_id'], 
                          "错误原因" : rxn['rules'].keys(),
                          "平衡情况" : rxns.check_mass_balance(),
                          "反应修正情况" : bounds_modification}
                infinite_info.append(infiniteInfo)
        model_control_info['nets']['修正反应情况'] = infinite_info


    def get_associated_substances(self, modelInfo, model, model_control_info, infinite_rxn, all_net):

        temp_net = []
        
        for netId in all_net:
            with model:
                object_rxn_id = add_demand(modelInfo, model, netId)
                model.objective = object_rxn_id
                if model.optimize().fluxes[object_rxn_id] <= 1e-6:
                    temp_net.append(netId)
        # nets_dict[','.join(infinite_rxn)] = temp_net  # 每个净物质处理后都会影响部分其他净物质，nets_dict字典记录了相关联的物质{'for_c': ['for_c', 'hco3_c', 'mmcoa__S_c', 'succoa_c', 'cbp_c', 'urea_c', 'allphn_c', 'agps23_c'], 'malcoa_c': ['malcoa_c', 'prpncoa_c', '3hpcoa_c', 'ppcoa_c', '3opcoa_c'], 'ag160_c': ['ag160_c', '2agpe160_c']}
        model_control_info["nets"]["净物质之间相关联情况"][f"{','.join(infinite_rxn)}"] = temp_net
        return temp_net
        # for k in temp_net:
        #     all_net.remove(k)
      

    def boundary_restoration(self, model_info, model):
        """
        Response to changes each cycle, restoring bounds in time
        """
        for rxn in model_info['reactions']:
            if rxn['net_modify'] == "true":
                model.reactions.get_by_id(rxn['id']).bounds = rxn['bounds']

    def net_control(self, model_info, model, check_model, model_control_info):
        """
        nets correction process
        """
        model_control_info["nets"], all_net = {}, []
        if not is_bio_exp_correct(model_info):
            model_control_info["biomass_info"] = "无法找到biomass方程,请检查后再试!!!"
            return 0
        rules = Rules()   
        atpm_id = self.get_atpm_rxn(model_info, model)
        set_model_objective(model_info, model, atpm_id)
        model_control_info["initial_rxn"] = model_info["initial_rxn"]
        model_control_info["initial_atpm_rxn"] = model_info["now_obj_rxn"]
        self.set_no_c_source_supply(model_info, model)
        if model.optimize().objective_value <= 1e-6:
            model_control_info["atps"]["净物质分数:"] = "100% (经过交换反应的质子设置,能量已经修复)"
        self.find_net_generation(model_info, model, model_control_info)
        all_net.extend(model_control_info["nets"]["净物质:"])
        model_control_info["nets"]["净物质之间相关联情况"] = {}
        for netId in model_control_info["nets"]["净物质:"]:
            with model:
                object_rxn_id = add_demand(model_info, model, netId)
                model.objective = object_rxn_id
                self.net_penalty_points(model_info, model, rules)
                grate_delete_rxn = self.close_max_value_net(model_info, model_control_info, model, object_rxn_id, netId)
                if len(grate_delete_rxn) == 0 :
                    continue
                infinite_rxn = self.add_back_reactions_net(model_info, model, check_model, grate_delete_rxn)
                self.get_net_infinite_info(model_info, model, model_control_info)
                temp_net = self.get_associated_substances(model_info, model, model_control_info, infinite_rxn, all_net)
                for k in temp_net:
                    all_net.remove(k)
            self.boundary_restoration(model_info, model)
            
