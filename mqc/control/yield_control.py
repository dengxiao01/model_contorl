# -*- coding: utf-8 -*-
import re 

from mqc.utils import *

class Yields():
    """"""

    def __init__(self) -> None:
        pass

    
    def is_autotrophic(self, model_info):
        """"""
        for rxn in model_info["reactions"]:
            if 'autotrophic' in rxn.keys(): 
                return True
        return False



    def check_yield_formula(self, model, metId):
        """"""
        calculate_met = ['C','H','O','N','P','S']
        met = model.metabolites.get_by_id(metId)
        if not set(met.elements.keys()).issubset(set(calculate_met)):
            return 0
        if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula or met.formula == 'null':
            return 0
        return 1


    def find_yield_generation(self, model_info, model, model_control_info):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux = {}
        for met in model_info['metabolites']:
            if re.search("C[A-Z]",str(met['formula'])) or re.search("C[\d]",str(met['formula'])): # Carbon followed by capital letters or numbers
                c_num += 1
                if self.check_yield_formula(model, met['id']) == 0:
                    print(met,': formula不明确,不计算')
                    continue
                
                with model:
                    objectiveId = add_demand(model_info, model, met['id'])
                    model.objective = objectiveId
                    write_flux_file(model_info, model, objectiveId)
                    max_rate = max_reduced_degree_rate(model, met['id'])
            
                    if max_rate == 0 or 'R' in met['formula'] or 'X' in met['formula']: # 含碳物质，最大得率如果是0不用算；formula含有[CHONPS]以外的不算
                        continue
                    print('met: ',met['id'],' object: ',model.slim_optimize(),'  max_rate: ',max_rate)
                    if round(model.slim_optimize(),2) > round(max_rate,2):
                        yield_generation.append(met['id'])
                        initial_yield_flux[f"{met['id']}"] = round(model.slim_optimize(),2)
        print('yield_generation:',yield_generation)
        if len(yield_generation) != 0:
            model_control_info["yields"]["score:"] = 0
        model_control_info["yields"]["含碳物质总数:"] = c_num
        model_control_info["yields"]["得率超过理论最大值数量:"] = len(yield_generation)
        model_control_info["yields"]["得率净物质:"] = yield_generation
        model_control_info["initial_yield_flux"] = initial_yield_flux





    def close_max_value_yield(self, model_info, model_control_info, model, object_rxn_id, yieldId, max_rate):
        """"""
        grate_delete_rxn, count = [], 0
        if round(model.slim_optimize(),2) <= round(max_rate,2):
            print(model,yieldId,'________得率_yes')
        while round(model.slim_optimize(),2) > round(max_rate,2): 
            count += 1
            print(count,' -- ',object_rxn_id)
            if count > 100:
                print('error: 得率无限循环')
                model_control_info["yields"] = "error: 得率无限循环"
                exit()
            write_flux_file(model_info, model, object_rxn_id)
            fba_rxn_dic = {}
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            for ids,v in need_fluxes.items():
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        if ids == object_rxn_id: 
                            continue
                        fba_rxn_dic[ids] = rxn['net_penalty_points']
            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()): #关掉每次结果中的所有最大值
                    if fba_rxn_dic[ids] != 0:  
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(model.reactions.get_by_id(ids),fba_rxn_dic[ids])
            max_rate = max_reduced_degree_rate(model, yieldId)
            print('close :  objective_value:',model.slim_optimize(),' max_rate: ',max_rate)
        return grate_delete_rxn

    def get_associated_substances(self, modelInfo, model, model_control_info, infinite_rxn, all_yield):
        """"""
        temp_yield = []
        
        for yieldId in all_yield:
            with model:
                object_rxn_id = add_demand(modelInfo, model, yieldId)
                model.objective = object_rxn_id
                temp_max_rate = max_reduced_degree_rate(model, yieldId)
                if round(model.slim_optimize(),2) <= round(temp_max_rate,2):
                    temp_yield.append(yieldId)
        # nets_dict[','.join(infinite_rxn)] = temp_net  # 每个净物质处理后都会影响部分其他净物质，nets_dict字典记录了相关联的物质{'for_c': ['for_c', 'hco3_c', 'mmcoa__S_c', 'succoa_c', 'cbp_c', 'urea_c', 'allphn_c', 'agps23_c'], 'malcoa_c': ['malcoa_c', 'prpncoa_c', '3hpcoa_c', 'ppcoa_c', '3opcoa_c'], 'ag160_c': ['ag160_c', '2agpe160_c']}
        model_control_info["yields"]["含碳物质依赖情况"][f"{','.join(infinite_rxn)}"] = temp_yield
        return temp_yield


    



    def yield_control(self, model_info, model, check_model, model_control_info):
        """"""
        all_yield = []   
        model_control_info["yields"]["score:"] = 1
        if self.is_autotrophic(model_info):
            return 0
        else:
            set_c_source_supply(model_info, model, 'yields')
            self.find_yield_generation(model_info, model, model_control_info)
            all_yield.extend(model_control_info["yields"]["得率净物质:"])
            model_control_info["yields"]["含碳物质依赖情况"] = {}
            model_control_info['yields']['模型修正情况'] = []
            for yieldId in model_control_info["yields"]["得率净物质:"]:
                print('yieldId:',yieldId)
                
                with model:
                    object_rxn_id = add_demand(model_info, model, yieldId)
                    model.objective = object_rxn_id 
                    
                    max_rate = max_reduced_degree_rate(model, yieldId)
                    print('objective_value: ',model.slim_optimize(),' max_rate: ',max_rate)
                    grate_delete_rxn = self.close_max_value_yield(model_info, model_control_info, model, object_rxn_id, yieldId, max_rate)
                    if len(grate_delete_rxn) == 0 :
                        continue
                    infinite_rxn = add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, 'yields', max_rate)
                    get_net_infinite_info(model_info, model, model_control_info, 'yields')
                    temp_net = self.get_associated_substances(model_info, model, model_control_info, infinite_rxn, all_yield)
                    for k in temp_net:
                        all_yield.remove(k)
                boundary_restoration(model_info, model, 'yields')
            model_control_info["final_yield_flux"] = get_final_yield_fluxes(model_info, model, model_control_info)
