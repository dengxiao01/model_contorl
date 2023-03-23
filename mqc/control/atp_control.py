# -*- coding: utf-8 -*-

from mqc.utils import *
from mqc.control.rules_control import Rules


class Atps():

    def __init__(self):
        """"""
        


    def modify_nadh(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        for ids,v in need_fluxes.items(): 
            for rxn in model_info['reactions']:
                if ids == rxn['id']:
                    rxns = model.reactions.get_by_id(ids)
                    if len(set(NADH) & set(rxn['reactants_mets'])) != 0 and v < 0 :
                        rxns.bounds = (0,1000)
                        rxn["bounds"] = [0,1000]
                        rxn["rules"]["modify_nadh_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                    if len(set(NADH) & set(rxn['products_mets'])) != 0 and v > 0 :
                        rxns.bounds = (-1000,0)
                        rxn["bounds"] = [-1000,0]
                        rxn["rules"]["modify_nadh_rxn"] = "true"
                        grate_delete_rxn.append(ids)


    def modify_pmf(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        for ids,v in need_fluxes.items(): 
            for rxn in model_info['reactions']:
                if ids == rxn['id']:
                    rxns = model.reactions.get_by_id(ids)
                    reactants_mets = [m.id for m in rxns.reactants]
                    products_mets = [m.id for m in rxns.products]
                    h_close = model_info["h_close"]
                    print(h_close)
                    if len(set(h_close) & set(reactants_mets)) != 0 and v < 0:
                        rxns.bounds = (0,1000)
                        rxn["bounds"] = [0,1000]
                        rxn["rules"]["modify_pmf_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                    if len(set(h_close) & set(products_mets)) != 0 and v > 0:
                        rxns.bounds = (-1000,0)
                        rxn["bounds"] = [-1000,0]
                        rxn["rules"]["modify_pmf_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                    # if rxn['id'] in model_info['exchange_rxns'] or rxn['id'] in model_info['transport_rxns']:
                    #     if (len(set(H_E) & set(reactants_mets)) != 0 or len(set(H_E) & set(rxn['reactants_mets'])) != 0) and v > 0:
                    #         rxns.bounds = (0,1000)
                    #         rxn["bounds"] = [0,1000]
                    #         rxn["rules"]["modify_pmf_rxn"] = "true"
                    #         grate_delete_rxn.append(ids)
                    #     if (len(set(H_E) & set(products_mets)) != 0 or len(set(H_E) & set(rxn['reactants_mets'])) != 0) and v < 0:
                    #         rxns.bounds = (-1000,0)
                    #         rxn["bounds"] = [-1000,0]
                    #         rxn["rules"]["modify_pmf_rxn"] = "true"
                    #         grate_delete_rxn.append(ids)


    def modify_atp(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        for ids,v in need_fluxes.items(): 
            for rxn in model_info['reactions']:
                if ids == rxn['id']:
                    if len(set(ATP_NAME) & set(rxn['reactants_mets'])) != 0 and len(set(ATP_NAME) & set(rxn['products_mets'])) != 0:
                        break 
                    rxns = model.reactions.get_by_id(ids)
                    if len(set(ATP_NAME) & set(rxn['reactants_mets'])) != 0 and v < 0 and 'right_chain_rxn' not in rxn['rules'].keys():
                        rxns.bounds = (0,1000)
                        rxn["bounds"] = [0,1000]
                        rxn["rules"]["modify_atp_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                    if len(set(ATP_NAME) & set(rxn['products_mets'])) != 0 and v > 0 and 'right_chain_rxn' not in rxn['rules'].keys():
                        rxns.bounds = (-1000,0)
                        rxn["bounds"] = [-1000,0]
                        rxn["rules"]["modify_atp_rxn"] = "true"
                        grate_delete_rxn.append(ids)



    def close_max_value_atp(self, model_info, model_control_info, model):
        """
        Turn off the maximum value in each result
        """
        grate_delete_rxn, count = [], 0
        if model.slim_optimize() <= 1e-6:
            print(model,'________ATP_yes')
        while model.slim_optimize() > 1e-6:
            count = count+1
            print(count)
            if count > 200:
                print('error: 能量无限循环,,模型中可能有额外底物输入或者错误反应')
                model_control_info["atps"] = "error: 能量无限循环,,模型中可能有额外底物输入或者错误反应"
                exit()
            atpm_id = model_control_info["initial_atpm_rxn"]["rxn_id"]
            write_flux_file(model_info, model, atpm_id)
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            if count > 50 : self.modify_nadh(model_info, model, need_fluxes, grate_delete_rxn)  
            if count > 100 : self.modify_pmf(model_info, model, need_fluxes, grate_delete_rxn)     
            if count > 150 : self.modify_atp(model_info, model, need_fluxes, grate_delete_rxn)
            fba_rxn_dic = {}
            for ids in need_fluxes.index:      # type(need_fluxes) : pandas.core.series.Series
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        fba_rxn_dic[ids] = rxn['net_penalty_points']
            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()):#关掉每次结果中的所有最大值
                    if ids != 'ADD_ATPM' and ids not in ATPM and fba_rxn_dic[ids]!= 0:  
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(model.reactions.get_by_id(ids),fba_rxn_dic[ids])
            print('close :',model.slim_optimize())
        return grate_delete_rxn
        

    def get_final_fluxes(self, model_info, model, check_model, model_control_info):
        """"""
        with model:
            for rxn in model_info['reactions']:
                if 'nadhs_modify' not in rxn.keys() and 'atps_modify' not in rxn.keys() and rxn['id'] not in ATPM:
                    rxns = model.reactions.get_by_id(rxn['id'])
                    if rxns.id != 'ADD_ATPM':
                        rxns.bounds = check_model.reactions.get_by_id(rxn['id']).bounds
            close_autotrophic_or_c_source(model_info, model, 'carbon_source')
            set_c_source_supply(model_info, model, 'atps')
            model_control_info["final_atp_flux"] = model.slim_optimize()
            atpm_id = model_control_info["initial_atpm_rxn"]["rxn_id"]
            write_flux_file(model_info, model, atpm_id)
            nadh_id = model_control_info["initial_nadh_rxn"]["rxn_id"]
            set_model_objective(model_info, model, nadh_id)
            model_control_info["final_nadh_flux"] = model.slim_optimize()
            # model_info['Carbon Source Information'].extend([f"Generated {model_control_info['final_nadh_flux']} NADH"])
            write_flux_file(model_info, model, nadh_id)
            if "ADD_ATPM" in model.reactions:
                model.reactions.remove('ADD_ATPM')
            model.reactions.remove('ADD_NADH')
            model.objective = model_info["initial_rxn"]["initial_rxn_id"]
            model_control_info["final_biomass_flux"] = model.slim_optimize()
            # model_control_info['Carbon Source Information'] = model_info['Carbon Source Information']
        # if "ADD_ATPM" in model.reactions:
        #     model.reactions.remove('ADD_ATPM')

  


    def atp_control(self, model_info, model, check_model, model_control_info):
        """"""
        model_control_info["atps"]["score:"] = 1
        num = 1e-6  
        atpm_id = get_atpm_rxn(model_info, model)
        set_model_objective(model_info, model, atpm_id)
        model_control_info["initial_atpm_rxn"] = model_info["now_obj_rxn"]
        if model.slim_optimize() <= 1e-6:
            model_control_info["atps"]["特殊修正:"] = "关闭碳源等供给源,能量得到修复"
        else:
            model_control_info["atps"]["score:"] = 0
        model_control_info['atps']['模型修正情况'] = []
        grate_delete_rxn = self.close_max_value_atp(model_info, model_control_info, model)
        add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, 'atps', num)
        get_net_infinite_info(model_info, model, model_control_info, 'atps')
        self.get_final_fluxes(model_info, model, check_model, model_control_info)