# -*- coding: utf-8 -*-

from mqc.utils import *
from mqc.control.rules_control import Rules


class Nadhs():

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


    def close_max_value_atp(self, model_info, model_control_info, model):
        """
        Turn off the maximum value in each result
        """
        grate_delete_rxn, count = [], 0
        if model.slim_optimize() <= 1e-6:
            print(model,'________NADH_yes')
        while model.slim_optimize() > 1e-6:
            count = count+1
            print(count)
            if count > 100:
                print('error: 还原力无限循环,,模型中可能有额外底物输入或者错误反应')
                model_control_info["nadhs"] = "error: 还原力无限循环,,模型中可能有额外底物输入或者错误反应"
                exit()
            atpm_id = model_info['now_obj_rxn']['now_obj_rxn_id']
            write_flux_file(model_info, model, atpm_id)
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            if count > 50 : self.modify_nadh(model_info, model, need_fluxes, grate_delete_rxn)  
            fba_rxn_dic = {}
            for ids in need_fluxes.index:      # type(need_fluxes) : pandas.core.series.Series
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        fba_rxn_dic[ids] = rxn['net_penalty_points']
            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()):#关掉每次结果中的所有最大值
                    if ids != 'ADD_NADH' and fba_rxn_dic[ids]!= 0:  
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(model.reactions.get_by_id(ids),fba_rxn_dic[ids])
            print('close :',model.slim_optimize())
        return grate_delete_rxn


    def get_final_fluxes(self, model_info, model, check_model, model_control_info):
        """"""
        with model:
            for rxn in model_info['reactions']:
                if 'nadhs_modify' not in rxn.keys():
                    rxns = model.reactions.get_by_id(rxn['id'])
                    if rxns.id != 'ADD_NADH':
                        rxns.bounds = check_model.reactions.get_by_id(rxn['id']).bounds
            close_autotrophic_or_c_source(model_info, model, 'carbon_source')
            set_c_source_supply(model_info, model, 'nadhs')
            model_control_info["final_nadh_flux"] = model.slim_optimize()
            print(model.slim_optimize())
            print(model.optimize().objective_value)
            nadh_id = model_info['now_obj_rxn']['now_obj_rxn_id']
            write_flux_file(model_info, model, nadh_id)
            model.reactions.remove('ADD_NADH')
            model.objective = model_info["initial_rxn"]["initial_rxn_id"]
            model_control_info["final_biomass_flux"] = model.slim_optimize()
        # model.reactions.remove('ADD_NADH')


    def nadh_control(self, model_info, model, check_model, model_control_info):
        """"""
        model_control_info["nadhs"]["score:"] = 1
        if not is_bio_exp_correct(model_info):
            model_control_info["biomass_info"] = "无法找到biomass方程,请检查后再试!!!"
            model_control_info["nadhs"]["score:"] = 0
            return 0
        set_model_boundary(model_info, model)
        model_control_info['boundary Information'] = model_info['boundary Information']
        num = 1e-6
        rules = Rules()   
        nadh_id = add_nadh_rxn(model)
        set_model_objective(model_info, model, nadh_id)
        model_control_info["initial_rxn"] = model_info["initial_rxn"]
        model_control_info["initial_nadh_rxn"] = model_info["now_obj_rxn"]
        close_autotrophic_or_c_source(model_info, model, 'carbon_source')
        close_autotrophic_or_c_source(model_info, model, 'autotrophic')
        if model.slim_optimize() <= 1e-6:
            model_control_info["nadhs"]["特殊修正:"] = "关闭碳源等供给源,还原力得到修复"
        else:
            model_control_info["nadhs"]["score:"] = 0
        rules.get_all_rules(model_info, model)
        check_rxn_balance(model_info, model)
        check_C_balance(model_info, model)
        net_penalty_points(model_info)
        model_control_info['nadhs']['模型修正情况'] = []
        grate_delete_rxn = self.close_max_value_atp(model_info, model_control_info, model)
        add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, 'nadhs', num)
        get_net_infinite_info(model_info, model, model_control_info, 'nadhs')
        # self.get_final_fluxes(model_info, model, check_model, model_control_info)