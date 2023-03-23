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



    def find_net_generation(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        r_met, net_generation, only_pro_mets, met_num = [], [], [], 0
        initial_net_flux = {}
        for rxn in model_info["reactions"]:  # Find the metabolites of the reactions with the lowest bound less than 0 in the exchange reaction
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"]
            if len(rxn['reactants_mets']) == 0 and len(rxn['products_mets']) != 0:
                only_pro_mets += rxn['products_mets']
        for met in model_info['metabolites']:
            met_num += 1
            if met['name'] not in (FREE_METS + only_pro_mets):
                if re.search("C[A-Z]",str(met['formula'])) or re.search("C[\d]",str(met['formula'])): # Carbon followed by capital letters or numbers         
                    with model:
                        objectiveId = add_demand(model_info, model, met['id'])
                        model.objective = objectiveId
                        if model.slim_optimize() > 1e-5 and met['name'] not in r_met:
                            net_generation.append(met['id'])
                            initial_net_flux[f"{met['id']}"] = model.slim_optimize()
        if len(net_generation) != 0:
            model_control_info["nets"]["score:"] = 0
        model_control_info["nets"]["代谢物总数:"] = met_num
        model_control_info["nets"]["净物质总数:"] = len(net_generation)
        model_control_info["nets"]["净物质:"] = net_generation
        model_control_info["initial_net_flux"] = initial_net_flux


    


    def close_max_value_net(self, model_info, model_control_info, model, object_rxn_id, netId):
        """
        Turn off the maximum value in each result
        """
        grate_delete_rxn = []
        count = 0
        if model.slim_optimize() <= 1e-5:
            print(model,'.......',object_rxn_id,'..........yes')
        while model.slim_optimize() > 1e-5:  
            count += 1
            print(count,' -- ',object_rxn_id)
            if count > 100:
                print('error: 模型中可能有额外底物输入')
                model_control_info["nets"] = "error: 模型中可能有额外底物输入"
                exit()
            write_flux_file(model_info, model, object_rxn_id)
            fba_rxn_dic = {}
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            for ids,v in need_fluxes.items():      # i : 'ATPM' 'ADD_ATPM'
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        rxns = model.reactions.get_by_id(ids)
                        reactants_mets = [m.id for m in rxns.reactants]
                        products_mets = [m.id for m in rxns.products]
                        if rxn['balance'] == "false":  # 生成目标物质的反应不能关闭
                            if (v > 0 and netId in products_mets) or (v < 0 and netId in reactants_mets):
                                continue 
                        if ids == object_rxn_id: # 当前结果文件里的DM目标反应不能参与下面的关闭
                            continue
                        fba_rxn_dic[ids] = rxn['net_penalty_points']  # 'ATPM':3  'ADD_ATPM':0

            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()):#关掉每次结果中的所有最大值
                    if ids != 'ADD_ATPM' and ids not in ATPM and fba_rxn_dic[ids] != 0:  
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(model.reactions.get_by_id(ids),fba_rxn_dic[ids])
            print('close :',model.slim_optimize())
        return grate_delete_rxn


    def get_associated_substances(self, modelInfo, model, model_control_info, infinite_rxn, all_net):

        temp_net = []
        
        for netId in all_net:
            with model:
                object_rxn_id = add_demand(modelInfo, model, netId)
                model.objective = object_rxn_id
                if model.slim_optimize() <= 1e-6:
                    temp_net.append(netId)
        # nets_dict[','.join(infinite_rxn)] = temp_net  # 每个净物质处理后都会影响部分其他净物质，nets_dict字典记录了相关联的物质{'for_c': ['for_c', 'hco3_c', 'mmcoa__S_c', 'succoa_c', 'cbp_c', 'urea_c', 'allphn_c', 'agps23_c'], 'malcoa_c': ['malcoa_c', 'prpncoa_c', '3hpcoa_c', 'ppcoa_c', '3opcoa_c'], 'ag160_c': ['ag160_c', '2agpe160_c']}
        model_control_info["nets"]["净物质之间相关联情况"][f"{','.join(infinite_rxn)}"] = temp_net
        return temp_net
        # for k in temp_net:
        #     all_net.remove(k)
      

    

    def net_control(self, model_info, model, check_model, model_control_info):
        """
        nets correction process
        """
        all_net, num = [], 1e-6
        model_control_info["nets"]["score:"] = 1
        self.find_net_generation(model_info, model, model_control_info)
        all_net.extend(model_control_info["nets"]["净物质:"])
        model_control_info["nets"]["净物质之间相关联情况"] = {}
        model_control_info['nets']['模型修正情况'] = []
        
        for netId in model_control_info["nets"]["净物质:"]:
            with model:
                object_rxn_id = add_demand(model_info, model, netId)
                model.objective = object_rxn_id   
                grate_delete_rxn = self.close_max_value_net(model_info, model_control_info, model, object_rxn_id, netId)
                if len(grate_delete_rxn) == 0 :
                    continue
                infinite_rxn = add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, 'nets', num)
                get_net_infinite_info(model_info, model, model_control_info, 'nets')
                temp_net = self.get_associated_substances(model_info, model, model_control_info, infinite_rxn, all_net)
                for k in temp_net:
                    all_net.remove(k)
            boundary_restoration(model_info, model, 'nets')
        model_control_info["final_net_flux"] = get_final_net_fluxes(model_info, model, model_control_info)
            
