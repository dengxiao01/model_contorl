# -*- coding: utf-8 -*-

import cobra


from mqc.control.model_control import ModelPreprocess
from mqc.control.preprocessing_control import Preprocess
from mqc.utils import *

class Biomasses():


    def __init__(self):
        """
        Get model information and check_model object.
        """

    def recover_bio_objective(self, model_info, model):
        """"""
        model.objective = model_info["initial_rxn"]["initial_rxn_id"]


    def test(self, model, check_model):
        """"""
        # model.reactions.get_by_id('EX_h_e').bounds = (0,1000)
        # model.reactions.get_by_id('FDR').bounds =(0,0)
        # model.reactions.get_by_id('G3POACTF').bounds = (0,0)
        # model.reactions.get_by_id('MADH').bounds = (0,0)
   

    def add_rely_biomass_rxn(self, model_info, model, check_model, model_control_info, initial_rxn_id):
        """"""
        temp, rely_rxn = [f"依赖于biomass的反应:"], []
        for rxn in model_info["reactions"]:
            if 'nadhs_modify' in rxn.keys() or 'atps_modify' in rxn.keys() or 'nets_modify' in rxn.keys() or 'yields_modify' in rxn.keys():
                with check_model:
                    check_model.objective = initial_rxn_id
                    check_model.reactions.get_by_id(rxn['id']).bounds = (0,0)
                    if check_model.slim_optimize() <= 1e-6:
                        rely_rxn.append(rxn['id'])
                        rxn["rely_biomass"] = "true"
            if 'nadhs_modify' not in rxn.keys() and 'atps_modify' not in rxn.keys() and 'nets_modify' not in rxn.keys() and 'yields_modify' not in rxn.keys() and 'carbon_source' not in rxn.keys():
                if rxn['id'] != initial_rxn_id:
                    model.reactions.get_by_id(rxn['id']).bounds = check_model.reactions.get_by_id(rxn['id']).bounds
        for rxnId in rely_rxn:
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(rxnId).bounds
        if rely_rxn:
            temp.extend(rely_rxn)
            write_biomass(model_control_info, temp, "预处理")


    def get_general_library(self, model_info):
        """"""
        if model_info['model_identifier'] == 'modelseed':
            general_library = cobra.io.read_sbml_model("mqc/summary/seed_1.xml")
        elif model_info['model_identifier'] == 'bigg':
            general_library = cobra.io.load_json_model("mqc/summary/MSGEM_end_with_unbalance_20220920.json")
        else:
            general_library = ''
        return general_library


    def gapfilling(self, model, demand_id, general_library):
        """"""
        with model :
            for rxns in model.reactions:  # 给系数扩大1000倍,添加的demand反应不改
                temp_dict = {}
                for met in rxns.metabolites:
                    temp_dict[met] = rxns.get_coefficient(met) * (1000)
                rxns.add_metabolites(temp_dict) # 相同时会覆盖掉原代谢物。乘以1000是为了pfba计算
            for rxns in general_library.reactions:
                if rxns not in model.reactions and 'BIOMASS' not in rxns.id and len(rxns.products) != 0:
                    model.add_reaction(rxns)  
            # model.reactions.get_by_id('EX_glc__D_e').bounds = (-1000,0)  #葡萄糖 （-1000，0）
            model.reactions.get_by_id(demand_id).bounds = (0,1)  #dna （0，1） 
            
            pfba_solution = pfba(model)  # 将模型和总库的反应融合到一起后，再次计算当前目标下的通量结果文件；随后把生成目标的途径添加到模型中
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            new_way_rxn = list(need_fluxes.keys())
        return new_way_rxn


    def add_gapfilling_rxn(self, model, model_info, demand_id, general_library):
        new_way_rxn = self.gapfilling(model, demand_id, general_library)
        num = 0
        if not new_way_rxn:
            gap_info = ["Not able to gapfill the response!!!"]
        else:
            need_add_rxn = []
            gap_info = ["gapfilling added reaction:"]
            for rxnId in new_way_rxn:
                if rxnId not in model.reactions and rxnId in general_library.reactions:
                    num += 1
                    rxns = general_library.reactions.get_by_id(rxnId)
                    model.add_reaction(rxns)
                    need_add_rxn.append(rxnId)
                    rxn_exp = [f"{rxns}"]
                    if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                        rxn_exp = [f"{rxns.id}: {rxns.build_reaction_string(use_metabolite_names=True)}"]
                    gap_info.extend(rxn_exp)
            if num == 0:
                gap_info.extend(["no gap response,since it already exists!"])
        return need_add_rxn, gap_info


    def add_gap_rxn(self, model, need_add_rxn, general_library):
        """"""
        for rxnId in need_add_rxn:
            rxns = general_library.reactions.get_by_id(rxnId)
            model.add_reaction(rxns)
      

    def check_bio_component_is_zero(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        no_synthesized_mets, num = [], 0
        for met in model.reactions.get_by_id(initial_rxn_id).reactants:
            num += 1
            with model:
                objectiveId = add_demand(model_info, model, str(met.id))
                model.objective = objectiveId
                ini_flux = model.slim_optimize()
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_synthesized_mets.append(str(met.name))
        temp = [f"Biomass consists of {num} components, {','.join(no_synthesized_mets)} is the rate of 0"]
        write_biomass(model_control_info, temp, "预处理")


    def gap_initial_zero_bio(self, model_info, model, model_control_info, initial_rxn_id, general_library):
        """"""
        no_synthesized_mets, num = [], 0
        temp = [f"Biomass initial rate is 0, check precursor components"]
        for met in model.reactions.get_by_id(initial_rxn_id).reactants:
            num += 1
            with model:
                objectiveId = add_demand(model_info, model, str(met.id))
                model.objective = objectiveId
                ini_flux = model.slim_optimize()
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    temp.extend([f"The {str(met.id)} rate of the biomass component material is 0, and try to gap:"])
                    no_synthesized_mets.append(str(met.id))
                    need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, objectiveId, general_library) 
                    temp.extend(gap_info)    
                    self.add_gap_rxn(model, need_add_rxn, general_library)
        write_biomass(model_control_info, temp, "预处理")
        return need_add_rxn, no_synthesized_mets


    

    def preprocess_initial_zero_bio(self, model_info, model, model_control_info, initial_rxn_id, general_library):
        """"""
        if model_info["initial_rxn"]["initial_rxn_flux"] == 0:
            model_control_info['biomasses']["score"] = 0
            need_add_rxn, no_synthesized_mets = self.gap_initial_zero_bio(model_info, model, model_control_info, initial_rxn_id, general_library)
            if need_add_rxn : self.add_gap_rxn(model, need_add_rxn, general_library)
            if model.slim_optimize() == 0:
                model_control_info['biomasses'] = "该模型无法模拟生长,所有小分子和大分子，都不能合成，请检查模型边界条件和单体合成途径"
                return -1
            else:
                model_control_info['biomasses'] = f"该模型生长速率为0,无法模拟生长,经检查后，前体物质{','.join(no_synthesized_mets)}不能合成,经过gap,可尝试添加{','.join(need_add_rxn)}反应来解决"
                model_info["initial_rxn"]["initial_rxn_flux"] = model.slim_optimize()
        return 0
     


    def get_adp_h_pi(self, model):
        """"""
        adp_identifier, h_identifier, pi_identifier = '', '', ''
        for met in model.metabolites:
            if met.id in ADP:
                adp_identifier = met.id 
            if met.id in H:
                h_identifier = met.id 
            if met.id in PI:
                pi_identifier = met.id 
        return adp_identifier, h_identifier, pi_identifier


    def is_missing_sustaining_energy(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        initial_rxns = model.reactions.get_by_id(initial_rxn_id)
        bio_text = model_info["initial_rxn"]["initial_rxn_exp"]
        missing_identifier, sustain, is_missing = [], [], 0
        products_mets = [met.id for met in model.reactions.get_by_id(initial_rxn_id).products]
        adp_identifier, h_identifier, pi_identifier = self.get_adp_h_pi(model)
        if len(set(ADP+H+PI) & set(products_mets)) == 0:
            sustain = ["biomass中生长相关维持(h、adp、pi)均缺失,请补足"]
            is_missing += 1
        elif len(set(ADP) & set(products_mets)) == 0 or len(set(H) & set(products_mets)) == 0 or len(set(PI) & set(products_mets)) == 0:
            if len(set(ADP) & set(products_mets)) == 0:  
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {adp_identifier}'
                missing_identifier.append(adp_identifier)
            elif len(set(H) & set(products_mets)) == 0:
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                missing_identifier.append(h_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {h_identifier}'
            elif len(set(PI) & set(products_mets)) == 0:
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                missing_identifier.append(pi_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {pi_identifier}'
            sustain = [f"biomass中生长相关维持能错误,有缺失,方程产物中添加{','.join(missing_identifier)}"]
            is_missing += 1
        if len(sustain) != 0:
            write_biomass(model_control_info, sustain, "预处理")
        initial_bio_value = model.slim_optimize()
        initial_biomass_rxn = initial_rxns.build_reaction_string(use_metabolite_names=False)
        if is_missing != 0:
            model_control_info['biomasses']["score"] = 0
            ini_bio_value = model_info["initial_rxn"]["initial_rxn_flux"]              
            sustain_result = [f"初始biomass生长速率(完善生长相关维持能前): {round(ini_bio_value,5)}"]
            sustain_result.extend([f"初始biomass方程(完善生长相关维持能前): {bio_text}"])
            sustain_result.extend([f"初始biomass生长速率(完善生长相关维持能后): {round(initial_bio_value,5)}"])
            sustain_result.extend([f"初始biomass方程(完善生长相关维持能后): {initial_biomass_rxn}"])
        else:
            sustain_result = [f"初始biomass生长速率(生长相关维持能完好): {round(initial_bio_value,5)}"]
            sustain_result.extend([f"初始biomass方程(生长相关维持能完好): {initial_biomass_rxn}"])
        return sustain_result


    def get_macromolecule(self, model, initial_rxn_id):
        """"""
        macromolecule = []
        for metId in model.reactions.get_by_id(initial_rxn_id).metabolites:
            met = model.metabolites.get_by_id(str(metId))
            if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula or met.formula == 'null':
                metName = met.name
                if any(metName.lower() in macro for macro in MACROMOLECULE) :  # 只要当前代谢物name在大分子列表里，就找到了大分子
                    macromolecule.append(str(metId))
        return macromolecule
    

    def get_smallmolecule(self, model_info, model, initial_rxn_id, macromolecule):
        """"""
        smallmolecule = []
        products_mets = [met.id for met in model.reactions.get_by_id(initial_rxn_id).products]
        for met in model_info['metabolites']:
            for metId in model.reactions.get_by_id(initial_rxn_id).metabolites:
                if met['id'] == metId.id:
                    if str(metId) not in macromolecule and str(metId) not in products_mets and met['name'] not in SMALLMOLECULE:
                        smallmolecule.append(str(metId))
        return smallmolecule


    def set_atpm_bounds(self, model_info, model):
        """"""
        for atpm_id in ATPM:
            if atpm_id in model_info['all_rxn_id']:
                model.reactions.get_by_id(atpm_id).bounds = (0,1000)
                for rxn in model_info['reactions']:
                    rxn['atpms_modify'] = "false"
                    if rxn['id'] == atpm_id:
                        rxn['atpms_modify'] = "true"


    def calculateResult(self, model, initial_rxn_id):
        """"""
        result, final = 0, 0
        pfba_solution = pfba(model)
        need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
        for ids,v in need_fluxes.items():
            v = round(v,5)
            rxns = model.reactions.get_by_id(ids)
            all_mets = [met.id for met in rxns.metabolites]
            if len(all_mets) == 1 and ids != initial_rxn_id and 'biomass' not in ids and 'BIOMASS' not in ids: #############################################用biomass列表替换
                mets = model.metabolites.get_by_id(all_mets[0])
                if pd.isna(mets.formula) or 'R' in mets.formula or 'X' in mets.formula or mets.formula == 'nan' or not mets.formula or mets.formula == 'null':
                    return all_mets[0]
                result += v * relative_molecular_mass(model, all_mets[0])
        try :
            final = abs(result/model.slim_optimize())
        except ZeroDivisionError :
            final = 0
        return final


    def get_initial_biomass(self, model_info, model, initial_rxn_id):
        """"""
        self.set_atpm_bounds(model_info, model)
        initial_biomass = self.calculateResult(model, initial_rxn_id)
        model_info['biomass']['initial_biomass'] = initial_biomass
        if type(initial_biomass) != str and (990 < initial_biomass < 1010):
            return 0



    def smallmolecule_iscoupling(self, model, model_info, small_demand_id, initial_rxn_id, metId, general_library, model_control_info):
        """"""
        iscoupling, close_coupling_rxn, need_add_rxn, temp_list = 0, [], [], []
        initial_flux = model.slim_optimize()
        print(initial_flux,type(initial_flux))
        if pd.isna(initial_flux) or initial_flux <= 1e-6 :
            temp_list = [f"{metId}速率为0,故进行gap:"] 
            need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, small_demand_id, general_library) 
            temp_list.extend(gap_info)
            write_biomass(model_control_info, temp_list, "小分子处理")
            if len(need_add_rxn) == 0:
                temp_list = [f"biomass生成中{metId}必须耦联大分子,gap不了,存在错误！！！请改正后再试"]
                write_biomass(model_control_info, temp_list, "小分子处理")
                return -1
            elif pd.isna(model.slim_optimize()):
                temp_list = [f"{metId}  gapfilling之后合成速率为nan,存在错误！！！请改正后再试"]
                write_biomass(model_control_info, temp_list, "小分子处理")
                return -1
            else:
                final_flux = model.slim_optimize()
                temp_list = [f"{metId}  gapfilling之后合成速率为: {round(final_flux,5)}"]
                write_biomass(model_control_info, temp_list, "小分子处理")
                return need_add_rxn, close_coupling_rxn
        while True :
            status = 0  # status状态符：只要本次循环没有发现biomass和大分子耦合反应，就退出循环
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            for r in list(need_fluxes.keys()): # 如果biomass反应也在结果中，就要解除大分子与biomass的耦联关系
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                if r == initial_rxn_id or (len(all_mets) == 1 and (pd.isna(met_formula) or 'R' in met_formula or 'X' in met_formula or met_formula == 'nan' or not met_formula or met_formula == 'null')):
                    if r != small_demand_id :  # 判断当前结果文件含有biomass耦合反应 或者 交换反应中大分子耦合反应
                        status +=  1
                        iscoupling += 1
                        model.reactions.get_by_id(r).bounds = (0,0)
                        close_coupling_rxn.append(r)
                        temp_list = [f"{small_demand_id}合成与{r}耦联,{metId}初始合成速率为： {round(initial_flux,5)},关闭{r}"]
                        write_biomass(model_control_info, temp_list, "小分子处理")
            if status == 0:
                break
        if iscoupling == 0 and initial_flux > 1e-6:
            return 0  
        else:
            final_flux2 = model.slim_optimize()
            if final_flux2 <= 1e-6 or pd.isna(final_flux2):
                temp_list = [f"小分子{metId}耦联的大分子全部关闭,生成速率为0, 进行gap:"]
                need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, small_demand_id, general_library) 
                temp_list.extend(gap_info)
                if len(need_add_rxn) == 0:
                    temp_list.extend([f"biomass生成中{metId}必须耦联大分子,gap不了,存在错误！！！请改正后再试"])
                    write_biomass(model_control_info, temp_list, "小分子处理")
                    return -1
                else :
                    final_flux = model.slim_optimize()
                    temp_list.extend([f"gapfilling之后合成速率为: {round(final_flux,5)}"])
                    write_biomass(model_control_info, temp_list, "小分子处理")
                    return need_add_rxn, close_coupling_rxn
            else:
                temp_list = [f"小分子{metId}耦联的大分子全部关闭,合成速率为{round(final_flux2,5)},不为0,则不进行gap"]
                write_biomass(model_control_info, temp_list, "小分子处理")
                return need_add_rxn, close_coupling_rxn


    def small_fix(self, model_info, model, initial_rxn_id, smallmolecule, general_library, model_control_info, check_model):
        """"""
        right_small_met, small_need_add_rxn, close_coupling_rxn = [], [], []
        for metId in smallmolecule:
            with model:
                small_demand_id = add_demand(model_info, model, str(metId))
                model.objective = small_demand_id
                if metId == 'cpd17042_c0':
                    write_flux_file(model_info, model, small_demand_id)
                small_info = self.smallmolecule_iscoupling(model, model_info, small_demand_id, initial_rxn_id, metId, general_library, model_control_info)
                if small_info == -1 :
                    return -1
                elif small_info == 0: 
                    right_small_met.append(metId)
                else:
                    small_need_add_rxn.extend(small_info[0])
                    close_coupling_rxn.extend(small_info[1])
            if close_coupling_rxn :
                model_info['biomass']['small_macro_modify'] = "yes"
                model_control_info['biomasses']["score"] = 0
                for r in close_coupling_rxn:
                    model.reactions.get_by_id(r).bounds = (0,0)
            if small_need_add_rxn:
                model_control_info['biomasses']["score"] = 0
                model_info['biomass']['small_macro_modify'] = "yes" 
                self.add_gap_rxn(model, small_need_add_rxn, general_library)
        temp = [f"以下小分子单体未耦联大分子，合成途径正确: {right_small_met}"]
        write_biomass(model_control_info, temp, "小分子处理")
        model.reactions.get_by_id(initial_rxn_id).bounds = check_model.reactions.get_by_id(initial_rxn_id).bounds


    def Normalized_macromolecule(self, model, model_info, metId, initialMacromoleculeResult, macro_demand_id, initial_rxn_id, normalized_macro_time, model_control_info, general_library):
        """"""
        close_coupling_rxn, initial_macromolecule_rxn, macromoleculeRxnId, need_add_rxn, close_coupling_rxn = [], '', '', [], []
        mets = model.metabolites.get_by_id(metId)
        while initialMacromoleculeResult :
            status = 0  # status状态符：只要本次循环没有发现biomass和大分子耦合反应，就退出循环
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            for r in list(need_fluxes.keys()): # 如果biomass反应也在结果中，就要解除大分子与biomass的耦联关系
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                if r == initial_rxn_id or (len(all_mets) == 1 and (pd.isna(met_formula) or 'R' in met_formula or 'X' in met_formula or met_formula == 'nan' or not met_formula or met_formula == 'null')):
                    if r != macro_demand_id :  # 判断当前结果文件含有biomass耦合反应 或者 交换反应中大分子耦合反应
                        status +=  1
                        model.reactions.get_by_id(r).bounds = (0,0)
                        close_coupling_rxn.append(r)
                        temp_Macromolecule = self.calculateResult(model, macro_demand_id) # 有耦合反应就把计算结果值赋给temp_Macromolecule
                        temp = [f"{macro_demand_id}合成与{r}耦联:{rxns}"]
                        temp.extend([f"关闭{r}后,{macro_demand_id}相对分子质量结果为{round(temp_Macromolecule,5)}"])
                        write_biomass(model_control_info, temp, "大分子处理")
                        if temp_Macromolecule == 0 :  # 判断当前结果是否为0，为0说明没有合成途径了，就gapfilling添加上途径
                            temp = [f"大分子全部关闭，{macro_demand_id}相对分子质量结果为0,进行gap:"]
                            need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, macro_demand_id, general_library) 
                            temp.extend(gap_info)
                            write_biomass(model_control_info, temp, "大分子处理")
                            if len(need_add_rxn) == 0:
                                temp = [f"biomass生成中{metId}必须耦联大分子,gap不了,存在错误,请改正"]
                                write_biomass(model_control_info, temp, "大分子处理")
                                return 0
                            break
                        else:  # 不为0时，可能是有途径，进入下一次循环，status为0退出循环；可能是还有其他的大分子反应，进入下一轮关
                            initialMacromoleculeResult = temp_Macromolecule
            if status == 0:
                break
        pfba_solution = pfba(model)
        need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6] 
        for r,v in need_fluxes.items():  # 获取合成大分子的反应以及代谢物系数
            rxns = model.reactions.get_by_id(r)
            products_mets = [m.id for m in rxns.products]
            reactants_mets = [m.id for m in rxns.reactants]
            coefficient_dict = {}
            if (mets.id in products_mets and v > 0) or (mets.id in reactants_mets and v < 0):
                macromoleculeRxnId = rxns.id
                initial_macromolecule_rxn = model.reactions.get_by_id(macromoleculeRxnId).build_reaction_string()
                for j in rxns.metabolites: # rxn是合成大分子的反应
                    # met_coefficient[str(j)] = rxn.get_coefficient(j)
                    if str(j) == mets.id: # 目标大分子系数不变
                        continue
                    coefficient_dict[j] = rxns.get_coefficient(j) * (1000/initialMacromoleculeResult) - rxns.get_coefficient(j) # 减去自己那一份，加回去就是所需的倍数
                rxns.add_metabolites(coefficient_dict)
                break 
        if normalized_macro_time > 1:
            temp = [f"第{normalized_macro_time}轮: {metId}方程: {initial_macromolecule_rxn}"]
            if status != 0 :
                temp.extend([f"第{normalized_macro_time}轮: {metId}耦联大分子{','.join(close_coupling_rxn)},且{','.join(close_coupling_rxn)}formula不明确!!!,无法计算初始{metId}的相对分子质量"])
            else :
                temp.extend([f"第{normalized_macro_time}轮: {mets.id}相对分子质量结果: {round(initialMacromoleculeResult,5)}"])
        else :
            temp = [f"初始{mets.id}方程: {initial_macromolecule_rxn}"]
            if status != 0 :
                temp.extend([f"{metId}耦联大分子{','.join(close_coupling_rxn)},且{','.join(close_coupling_rxn)}formula不明确!!!,无法计算初始{metId}的相对分子质量"])
            else :
                temp.extend([f"初始{metId}相对分子质量结果: {round(initialMacromoleculeResult,5)}"])
        write_biomass(model_control_info, temp, "大分子处理")
        return macromoleculeRxnId, need_add_rxn, close_coupling_rxn



    def recover_macromolecule_rxn(self, model, model_info, macromoleculeRxn, need_add_rxn, general_library, need_close_coupling_rxn, initial_rxn_id, model_control_info, check_model):
        """"""
        for k,v in macromoleculeRxn.items():
            model_control_info['biomasses']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            model.reactions.get_by_id(k).reaction = v
        if need_add_rxn:
            model_control_info['biomasses']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            self.add_gap_rxn(model, need_add_rxn, general_library)
        for r in need_close_coupling_rxn:
            model_control_info['biomasses']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            model.reactions.get_by_id(r).bounds = (0,0)
        model.reactions.get_by_id(initial_rxn_id).bounds = check_model.reactions.get_by_id(initial_rxn_id).bounds



    def macro_fix(self, model_info, model, initial_rxn_id, macromolecule, general_library, model_control_info, check_model):
        """"""
        gapfilling_rxn, need_close_coupling_rxn, macromoleculeRxn = [], [], {}
        for metId in macromolecule:
            normalized_macro_time, while_time = 0, 0
            with model:
                macro_demand_id = add_demand(model_info, model, str(metId))
                model.objective = macro_demand_id
                try :
                    model.reactions.get_by_id(macro_demand_id).bounds = check_model.reactions.get_by_id(macro_demand_id).bounds # 大分子反应边界还原
                except :
                    print(f'{model}没有{macro_demand_id}反应,不需要还原')
                print(f'{metId}_objective_value: ',model.slim_optimize())
                initialMacromoleculeResult = self.calculateResult(model, macro_demand_id)
                if type(initialMacromoleculeResult) != str and (990 < initialMacromoleculeResult < 1010) :
                    continue
                else:
                    normalized_macro_time += 1
                    macromoleculeInfo = self.Normalized_macromolecule(model, model_info, metId, initialMacromoleculeResult, macro_demand_id,initial_rxn_id,normalized_macro_time, model_control_info, general_library) 
                    if macromoleculeInfo == -1 :
                        return -1
                    elif macromoleculeInfo == 0:
                        continue
                    else :
                        macromoleculeRxnId, need_add_rxn, close_coupling_rxn = macromoleculeInfo
                    finalMacromoleculeResult = self.calculateResult(model, macro_demand_id)
                    gapfilling_rxn.extend(need_add_rxn)
                    need_close_coupling_rxn.extend(close_coupling_rxn)
                    while type(finalMacromoleculeResult) == str or finalMacromoleculeResult < 990 or finalMacromoleculeResult > 1010 :
                        normalized_macro_time += 1
                        macromoleculeInfo = self.Normalized_macromolecule(model, model_info, metId, finalMacromoleculeResult, macro_demand_id,initial_rxn_id,normalized_macro_time, model_control_info, general_library)
                        if macromoleculeInfo == -1 :
                            return -1
                        elif macromoleculeInfo == 0:  # 只是退出了第一个while
                            while_time += 1
                            break
                        else :
                            macromoleculeRxnId, need_add_rxn, close_coupling_rxn = macromoleculeInfo
                        finalMacromoleculeResult = self.calculateResult(model, macro_demand_id)
                        gapfilling_rxn.extend(need_add_rxn)
                        need_close_coupling_rxn.extend(close_coupling_rxn)
                    if while_time == 1 :
                        continue
                    finalRxn = model.reactions.get_by_id(macromoleculeRxnId).build_reaction_string(use_metabolite_names=False)
                    final_rxn_id = model.reactions.get_by_id(macromoleculeRxnId).build_reaction_string()
                    macromoleculeRxn[macromoleculeRxnId] = final_rxn_id
                    temp = [f"修正后{metId}方程: {finalRxn}"]
                    temp.extend([f"修正后{metId}相对分子质量结果: {round(finalMacromoleculeResult,5)}"])
                    write_biomass(model_control_info, temp, "大分子处理")
        self.recover_macromolecule_rxn(model, model_info, macromoleculeRxn, gapfilling_rxn, general_library, need_close_coupling_rxn, initial_rxn_id, model_control_info, check_model)



    def Normalized_biomass(self, model, initial_rxn_id, final_biomass):
        """"""
        coefficient_dict = {}
        if final_biomass == 0:
            return -1
        rxns = model.reactions.get_by_id(initial_rxn_id)
        for j in rxns.metabolites:
            if str(j) == 'biomass_c': # 如果方程有biomass_c，那么biomass_c系数不变##########################改biomass列表
                continue
            coefficient_dict[j] = rxns.get_coefficient(j) * (1000/final_biomass) - rxns.get_coefficient(j)
        rxns.add_metabolites(coefficient_dict)
        return 0


    def biomass_fix(self, model_info, model, initial_rxn_id, model_control_info):
        """"""
        middle_bio_value = model.slim_optimize()
        modify_nor_bio = final_biomass = self.calculateResult(model, initial_rxn_id)
        modify_nor_bio_eq = model.reactions.get_by_id(initial_rxn_id).build_reaction_string(use_metabolite_names=False)
        model_info['biomass']['middle_bio_value'] = middle_bio_value
        model_info['biomass']['modify_nor_bio_eq'] = modify_nor_bio_eq
        model_info['biomass']['modify_nor_bio'] = modify_nor_bio
        if type(final_biomass) != str and (990 < final_biomass < 1010):
            temp = [f"小分子大分子处理后biomass为1g"]
        else:
            temp = [f"小分子大分子处理后biomass不为1g"]
        write_biomass(model_control_info, temp, "biomass处理")
        if type(final_biomass) == str:
            temp = [f"最终biomass仍然耦联大分子{final_biomass}"]
            write_biomass(model_control_info, temp, "biomass处理")
            return -1
        while final_biomass < 990 or final_biomass > 1010:
            get_final = self.Normalized_biomass(model, initial_rxn_id, final_biomass)
            if get_final == -1:
                temp = [f"最终通量final为0,出错"]
                write_biomass(model_control_info, temp, "biomass处理")
                return -1
            final_biomass = self.calculateResult(model, initial_rxn_id)
        final_biomass_equation = model.reactions.get_by_id(initial_rxn_id).build_reaction_string(use_metabolite_names=False)
        final_biomass_equation_id = model.reactions.get_by_id(initial_rxn_id).build_reaction_string()
        model_info['biomass']['final_biomass'] = final_biomass
        model_info['biomass']['final_biomass_equation'] = final_biomass_equation
        model.reactions.get_by_id(initial_rxn_id).reaction = final_biomass_equation_id
                

    def recover_atpm_bounds(self, model_info, model, check_model):
        """"""
        for rxn in model_info['reactions']:
            if 'atpms_modify' in rxn.keys() and rxn['atpms_modify'] == "true":
                model.reactions.get_by_id(rxn['id']) == check_model.reactions.get_by_id(rxn['id'])

    def get_exchange(self, model, model_control_info):
        pfba_solution = pfba(model)
        need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
        temp = [f"模型计算条件:"]
        for r,v in need_fluxes.items():
            rxn = model.reactions.get_by_id(r)
            all_mets = [m.id for m in rxn.metabolites]
            # products_mets = [m.id for m in rxn.products]
            if len(all_mets) == 1:
                mets = model.metabolites.get_by_id(all_mets[0])
                temp.extend([f"{rxn}  {mets.formula}  {round(v,5)}  {rxn.bounds}"])
        write_biomass(model_control_info, temp, "模型计算条件")


    def get_result(self, model_info, model, sustain_result, model_control_info):
        """"""
        final_bio_value = model.slim_optimize()
        initial_biomass = model_info['biomass']['initial_biomass']
        if type(initial_biomass) == str:
            sustain_result.extend([f"biomass耦联大分子{initial_biomass}，且{initial_biomass}的formula不明确!!!,无法计算初始biomass的相对分子质量"])
        else:
            sustain_result.extend([f"初始biomass相对分子质量结果: {round(initial_biomass,5)}"])
        if 'small_macro_modify' in model_info['biomass'].keys():
            modify_nor_bio = model_info['biomass']['modify_nor_bio']
            modify_nor_bio_eq = model_info['biomass']['modify_nor_bio_eq']
            middle_bio_value = model_info['biomass']['middle_bio_value']
            if type(modify_nor_bio) == str:
                sustain_result.extend([f"修改小分子、大分子反应后biomass方程: {modify_nor_bio_eq}"])
                sustain_result.extend([f"biomass耦联大分子{modify_nor_bio}，且{modify_nor_bio}的formula不明确!!!,无法计算修改小分子、大分子反应后biomass的相对分子质量"])
            else:
                sustain_result.extend([f"修改小分子、大分子反应后biomass方程: {modify_nor_bio_eq}"])
                sustain_result.extend([f"修改小分子、大分子反应后biomass相对分子质量结果: {round(modify_nor_bio,5)}"])
            sustain_result.extend([f"修改小分子、大分子反应后biomass生长速率为: {round(middle_bio_value,5)}"])
        final_biomass_equation = model_info['biomass']['final_biomass_equation']
        final_biomass = model_info['biomass']['final_biomass'] 
        sustain_result.extend([f"修改biomass后,最终biomass方程: {final_biomass_equation}"])
        sustain_result.extend([f"修改biomass后,最终biomass相对分子质量结果: {round(final_biomass,5)}"])
        sustain_result.extend([f"最终biomass生长速率: {round(final_bio_value,5)}"])
        write_biomass(model_control_info, sustain_result, "biomass处理")


    def get_final_fluxes(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        atpm_id = get_atpm_rxn(model_info, model)
        set_model_objective(model_info, model, atpm_id)
        final_atp = model.slim_optimize()
        write_flux_file(model_info, model, atpm_id)
        if "ADD_ATPM" in model.reactions:
            model.reactions.remove('ADD_ATPM')
        nadh_id = add_nadh_rxn(model)
        set_model_objective(model_info, model, nadh_id)
        final_nadh = model.slim_optimize()
        write_flux_file(model_info, model, nadh_id)
        model.reactions.remove('ADD_NADH')
        final_netld_flux = get_final_net_fluxes(model_info, model, model_control_info)
        final_yield_flux = get_final_yield_fluxes(model_info, model, model_control_info)
        model.objective = initial_rxn_id
        final_biomass_flux = model.slim_optimize()
        model_control_info["final_flux"] = {'final_atp':final_atp,
                                            'final_nadh':final_nadh,
                                            'final_netld':final_netld_flux,
                                            'final_yield':final_yield_flux,
                                            'final_biomass':final_biomass_flux}


    def biomass_control(self, model_info, model, check_model, model_control_info):
        """"""
        print('biomass start')
        model_control_info['biomasses']["score"] = 1
        model_control_info['biomasses']["预处理:"] = []
        model_control_info['biomasses']["小分子处理:"] = []
        model_control_info['biomasses']["大分子处理:"] = []
        model_control_info['biomasses']["biomass处理:"] = []
        model_control_info['biomasses']["模型计算条件:"] = []
        model_info['biomass'] = {}
        self.recover_bio_objective(model_info, model)
        # self.test(model, check_model)
        initial_rxn_id = model_info["initial_rxn"]["initial_rxn_id"]
        self.add_rely_biomass_rxn(model_info, model, check_model, model_control_info, initial_rxn_id)
        general_library = self.get_general_library(model_info)     
        self.check_bio_component_is_zero(model_info, model, model_control_info, initial_rxn_id)
        if self.preprocess_initial_zero_bio(model_info, model, model_control_info, initial_rxn_id, general_library) == -1:
            return 0
        sustain_result = self.is_missing_sustaining_energy(model_info, model, model_control_info, initial_rxn_id)
        macromolecule = self.get_macromolecule(model, initial_rxn_id)
        smallmolecule = self.get_smallmolecule(model_info, model, initial_rxn_id, macromolecule)
        if self.get_initial_biomass(model_info, model, initial_rxn_id) == 0:
            model_control_info['biomasses']['biomass is ok?'] = 'yes'
            return 0
        if self.small_fix(model_info, model, initial_rxn_id, smallmolecule, general_library, model_control_info, check_model) == -1:
            return 0
        if self.macro_fix(model_info, model, initial_rxn_id, macromolecule, general_library, model_control_info, check_model) == -1:
            return 0
        if self.biomass_fix(model_info, model, initial_rxn_id, model_control_info) == -1:
            return 0
        self.recover_atpm_bounds(model_info, model, check_model)
        self.get_result(model_info, model, sustain_result, model_control_info)
        self.get_exchange(model, model_control_info)
        self.get_final_fluxes(model_info, model, model_control_info, initial_rxn_id)





        
