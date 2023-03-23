# -*- coding: utf-8 -*-

from cobra import Model, Reaction, Metabolite
import pandas as pd
import math 
import re 
from cobra.flux_analysis import pfba
from mqc.defaults import *

def get_model_file():
    """
    Get the model file.

    Return
    ------
    str
        model file

    """
    file = "mqc/test_data/bigg_data/iJN746.xml"
    return file


def add_rxn(model, rxn_name, rxn_expression):
    """
    Add the reaction to the model.

    Parameters
    ----------
    rxn_name : str
        Reactive name —— "ADD_ATPM"
    rxn_expression : str
        Reactive expression —— "atp_c + h2o_c  --> adp_c + h_c + pi_c"
    model : cobra.Model
 
    Return
    ------
    cobra.Model
        The model after adding the reaction

    """
    reaction = Reaction(rxn_name)
    model.add_reactions([reaction])
    reaction.build_reaction_from_string(rxn_expression)
    return model 


def add_atpm_rxn(model):
    """
    Add atpm reaction to model.

    Parameters
    ----------
    model : cobra.Model

    Return
    ------
    cobra.Model
        Model after adding atpm 

    """
    for met in model.metabolites:
        if met.id in ATP : atp_c = met.id
        if met.id in H : h_c = met.id 
        if met.id in ADP : adp_c = met.id 
        if met.id in H2O : h2o_c = met.id 
        if met.id in PI : pi_c = met.id 
    rxn_name = "ADD_ATPM"
    rxn_exp = "{} + {} --> {} + {} + {}".format(atp_c, h2o_c, adp_c, h_c, pi_c)
    add_rxn(model, rxn_name, rxn_exp)
    return rxn_name 


def add_nadh_rxn(model):
    """
    Add nadh reaction to model.

    Parameters
    ----------
    model : cobra.Model

    Return
    ------
    cobra.Model
        Model after adding nadh 

    """
    # for met in model.metabolites:
    #     if met.id in NADH_ID : nadh_c = met.id
    #     if met.id in H : h_c = met.id 
    #     if met.id in NAD_ID : nad_c = met.id 
    for met in model.metabolites:
        if met.id in NADH_ID : nadh_c = met.id
        if met.id in H : h_c = met.id 
        if met.id in NAD_ID : nad_c = met.id 
    rxn_name = "ADD_NADH"
    rxn_exp = "{} + {} --> {}".format(nadh_c, h_c, nad_c)
    add_rxn(model, rxn_name, rxn_exp)
    return rxn_name 




def get_atpm_rxn(model_info, model):
    """
    get atpm reactions
    """
    if len(set(model_info["all_rxn_id"]) & set(ATPM)) == 0:
        atpm_id = add_atpm_rxn(model)
    else:
        atpm_id = [i for i in ATPM if i in model_info["all_rxn_id"]][0]
    return atpm_id


def set_initial_obj(model):
    """
    Find the biomass equation and set it as the target.

    Parameters
    ----------
    model : cobra.Model

    Return
    ------
    str
        biomass equation ID

    """
    for rxn in model.reactions:
        all_mets = [m.id for m in rxn.metabolites]
        if len(all_mets) != 1 and ('bio' in rxn.id or 'biomass' in rxn.id or 'BIO' in rxn.id or 'BIOMASS' in rxn.id):
            return rxn.id
        else :
            return ''



def set_initial_obj2(model, rxnId):
    """
    Find the biomass equation and set it as the target.

    Parameters
    ----------
    rxnId : str
        Initial goals for model setup

    Notes
    -----
    The initial goal is biomass_c ⇌, Re-find the new biomass equation from the model

    """
    rxn = model.reactions.get_by_id(rxnId)
    all_mets = [m.id for m in rxn.metabolites]
    if len(all_mets) == 1 and ('bio' in rxn.id or 'biomass' in rxn.id or 'BIO' in rxn.id or 'BIOMASS' in rxn.id):
        for rxns in model.reactions:
            products_mets = [m.id for m in rxns.products]
            if all_mets[0] in products_mets:
                return rxns.id
            else:
                return ''
    else:
        return rxnId

def add_demand(modelInfo, model, metId):
    """
    Add demand reaction.

    Parameters
    ----------
    modelInfo : dict
        Information about the model
    metId : str
        The metabolite ID that needs to be added to the demand reaction

    Notes
    -----
    If there is an exchange reaction containing the metabolite ID in the model, set it directly as the target, otherwise add the demand reaction

    """
    mets = model.metabolites.get_by_id(metId) 
    flag, objectiveId = 0, ''
    for rxn in modelInfo["reactions"]:
        rxns = model.reactions.get_by_id(rxn['id'])
        all_mets = [met.id for met in rxns.metabolites]
        if ((rxn['reactants_mets'] and rxn["bounds"][1] > 0) or (rxn['products_mets'] and rxn['bounds'][0] < 0)) and metId in all_mets and len(all_mets) == 1:
            objectiveId = rxn["id"]
            flag += 1
    if flag == 0:
        if 'DM_' + metId in modelInfo["all_rxn_id"]:
            objectiveId = f"DM_{metId}"
        else:
            try:
                demand_rxn = model.add_boundary(mets, type = 'demand') # 添加demand反应
                objectiveId = demand_rxn.id
            except:
                objectiveId = 'DM_' + metId
    return objectiveId


def set_model_objective(modelInfo, model, rxnId):
    """
    Set the incoming response ID as the target and update modelInfo.
    """
    # now_obj_rxn_info = []
    model.objective = rxnId
    now_obj_rxn_flux = model.slim_optimize()
    if modelInfo['model_identifier'] == 'modelseed'or modelInfo['model_identifier'] == 'metaNetx':
        now_obj_rxn_exp = model.reactions.get_by_id(rxnId).build_reaction_string(use_metabolite_names=True)
    else:
        now_obj_rxn_exp = model.reactions.get_by_id(rxnId).build_reaction_string()
    NowObjRxnInfo = {f"rxn_id" : rxnId,
                    f"rxn_flux" : now_obj_rxn_flux,
                    f"rxn_exp" : now_obj_rxn_exp}
    # now_obj_rxn_info.append(NowObjRxnInfo)
    modelInfo["now_obj_rxn"] = NowObjRxnInfo


def is_bio_exp_correct(modelInfo):
    """
    Judging whether the biomass response of the model is correctly represented.
    """
    # id_index = modelInfo["initial_rxn"].index("initial_rxn_id")
    if not modelInfo["initial_rxn"]["initial_rxn_id"]:
        return 0
    return 1


def relative_molecular_mass(model, metid):
    """
    Judging whether the biomass response of the model is correctly represented.
    """
    met = model.metabolites.get_by_id(metid)
    mass = 0
    if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula:
        return 0
    for e,en in met.elements.items():
        if e == 'C': mass += 12*en       
        if e == 'H': mass += 1*en         
        if e == 'O': mass += 16*en        
        if e == 'N': mass += 14*en        
        if e == 'P': mass += 31*en        
        if e == 'S': mass += 32*en       
        if e == 'Ca': mass += 40*en         
        if e == 'Mn': mass += 55*en         
        if e == 'Fe': mass += 56*en     
        if e == 'Zn': mass += 65*en
        if e == 'Na': mass += 23*en
        if e == 'Mg': mass += 24*en
        if e == 'Al': mass += 27*en 
        if e == 'Si': mass += 28*en 
        if e == 'Cl': mass += 35.5*en 
        if e == 'K': mass += 39*en 
        if e == 'Cu': mass += 64*en 
        if e == 'Ag': mass += 108*en 
        if e == 'I': mass += 127*en 
        if e == 'Au': mass += 197*en 
        if e == 'Ba': mass += 137*en 
    return mass


def reduced_degree(model, metId):
    """Calculating the degree of reduction of a metabolite's single carbon"""
    met = model.metabolites.get_by_id(metId)
    degree = 0
    for e,en in met.elements.items():
        if e == 'C':
            degree = degree  + 4*en
        if e == 'H':
            degree = degree  + 1*en
        if e == 'O':
            degree = degree  + (-2)*en
        if e == 'N':
            degree = degree  + (-3)*en
        if e == 'P':
            degree = degree  + (5)*en
        if e == 'S':
            degree = degree  + (6)*en
    if pd.isna(met.charge) or met.charge == 'nan' or not met.charge:
        met.charge = 0
    degree = degree - met.charge
    return degree

def max_reduced_degree_rate(model, metId):
    """Calculate the theoretical yield of the maximum reduction degree of the product"""
    c_met_degree, max_rate = 0, 0
    pfba_solution = pfba(model)  
    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
    
    degree = reduced_degree(model, metId)
    if round(degree,5) != 0:
        for r,v in need_fluxes.items() :
            rxn = model.reactions.get_by_id(r)  
            all_mets = [m.id for m in rxn.metabolites]
            mets = model.metabolites.get_by_id(all_mets[0])
            if len(all_mets) == 1 and metId not in all_mets and 'C' in mets.elements.keys():
                if (rxn.reactants and rxn.lower_bound < 0) or (rxn.products and rxn.upper_bound > 0):
                    c_met_degree += (abs(round(v,5)) * reduced_degree(model, all_mets[0]))
        max_rate = c_met_degree/degree
    else:
        max_rate = 0
    return max_rate


def find_autotrophic_and_carbon_source(model_info, model):
    """"""
    for rxn in model_info["reactions"]:
        rxns = model.reactions.get_by_id(rxn['id'])
        if rxn["id"] in model_info["exchange_rxns"]:
            if (rxn['reactants_mets'] and rxn["bounds"][0] < 0) or (rxn['products_mets'] and rxn['bounds'][1] > 0):
                r_met_name = rxn["all_mets"][0] 
                r_met = list(rxns.metabolites.keys())[0]
                if any(r_met.id.startswith(source_name) for source_name in AUTOTROPHIC_SOURCE) or str(r_met.formula).find('Z') != -1 or r_met_name in PHOTON_NAME:
                    rxn["autotrophic"] = "true"
                if re.search("C[A-Z]",str(r_met.formula)) or re.search("C[\d]",str(r_met.formula)) or r_met_name in C_LIST or str(r_met.formula).find('R') != -1 or str(r_met.formula).find('X') != -1 or str(r_met.formula) == 'null': # C后边是大写字母或者数字
                    if r_met_name not in CO2_NAME or r_met_name not in HCO3_NAME:
                        rxn["carbon_source"] = "true"


def normalize_all_rxn_bounds(model_info, model):
    """change the reaction boundary outside the model exchange reaction to (-1000, 1000) within this range"""
    for rxn in model_info["reactions"]:
        rxns = model.reactions.get_by_id(rxn['id'])
        if rxn["bounds"][0] < 0 and rxn["bounds"][1] > 0:
            rxns.bounds = (-1000,1000)
            rxn["bounds"] = [-1000,1000]
        if rxn["bounds"][0] >= 0 and rxn["bounds"][1] > 0:
            rxns.bounds = (0,1000)
            rxn["bounds"] = [0,1000]
        if rxn["bounds"][0] < 0 and rxn["bounds"][1] <= 0:
            rxns.bounds = (-1000,0)
            rxn["bounds"] = [-1000,0]
            
    
                     

def close_autotrophic_or_c_source(model_info, model, flag):
    """Set no autotrophic or C source supply"""
    for rxn in model_info["reactions"]:
        if flag in rxn.keys():
            if rxn['reactants_mets'] and rxn["bounds"][0] < 0:
                model.reactions.get_by_id(rxn['id']).bounds = (0,1000)
                rxn["bounds"] = [0,1000]
            if rxn['products_mets'] and rxn['bounds'][1] > 0:
                model.reactions.get_by_id(rxn['id']).bounds = (-1000,0)
                rxn["bounds"] = [-1000,0]           



def set_model_boundary(model_info, model):
    """"""
    carbon_source_rxn, carbon_source_rxn2, num = [], [], 0
    model_info['boundary Information'], model_info['need_set_carbon_source'] = [], {}
    initial_rxn_id = model_info["initial_rxn"]["initial_rxn_id"]
    initial_bio_rxn = model.reactions.get_by_id(initial_rxn_id)
    if initial_bio_rxn.bounds != (0,1000):  
        model_info['boundary Information'].extend([f"The initial biomass boundary is {initial_bio_rxn.bounds}, changed to (0, 1000)"])
        initial_bio_rxn.bounds = (0,1000)
    else:
        model_info['boundary Information'].extend([f"The initial biomass boundary is (0, 1000),is ok"])
    normalize_all_rxn_bounds(model_info, model)
    find_autotrophic_and_carbon_source(model_info, model)
    for rxn in model_info["reactions"]:
        if 'carbon_source' in rxn.keys():
            carbon_source_rxn.append(rxn['id'])
            rxns = model.reactions.get_by_id(rxn['id'])   
            r_met = list(rxns.metabolites.keys())[0]
            if rxns.bounds != (-1000,1000) and str(r_met.formula) != 'null':
                num += 1
                break 
    if num == 0:
        model_info['boundary Information'].extend(["Wrong initial model bounds, all carbon source bounds are (-1000, 1000)"])
        glucose_rxnId = model_info['glucose_rxnId']
        if glucose_rxnId:
            model_info['need_set_carbon_source'][glucose_rxnId] = (-10,0)
            model_info['boundary Information'].extend([f"glucose present,Use glucose({glucose_rxnId}) as the carbon source,and supply 10"])
        else:
            model_info['need_set_carbon_source'][carbon_source_rxn[0]] = (-10,0)
            model_info['boundary Information'].extend([f"There is no glucose,Use model default carbon source {carbon_source_rxn[0]} as the carbon source,and supply 10"])
    else:  
        for rxn in model_info["reactions"]:
            rxns = model.reactions.get_by_id(rxn['id'])
            if "carbon_source" in rxn.keys():
                model_info['need_set_carbon_source'][rxn['id']] = rxns.bounds
                if rxns.lower_bound < -100:      
                    model_info['need_set_carbon_source'][rxn['id']] = (-10,0)
                    carbon_source_rxn2.append(rxn['id'])
        if carbon_source_rxn2:
            try:
                model_info['boundary Information'].extend([f"The initial boundary of the model is wrong, change {carbon_source_rxn2[0],carbon_source_rxn2[1],carbon_source_rxn2[2]} and other {len(carbon_source_rxn2)} carbon sources whose lowest boundary is less than -100 to (-10, 0), that is, 10 supply"])
            except:
                model_info['boundary Information'].extend([f"The initial boundary of the model is wrong, change {carbon_source_rxn2[0]} and other {len(carbon_source_rxn2)} carbon sources whose lowest boundary is less than -100 to (-10, 0), that is, 10 supply"])        
        else:
            model_info['boundary Information'].extend([f"The initial model boundary is correct, a total of {len(carbon_source_rxn)} carbon sources are supplied"])


def set_c_source_supply(model_info, model, flag):
    """"""
    for rxnId, rxn_bounds in model_info['need_set_carbon_source'].items():
        for rxn in model_info["reactions"]:
            if rxnId == rxn['id']:
                model.reactions.get_by_id(rxnId).bounds = rxn_bounds
                if flag == 'yields':
                    rxn['bounds'] = rxn_bounds
    


def rxn_no_annotation(model_info, model):
    for rxn in model_info['reactions']:
        rxns = model.reactions.get_by_id(rxn["id"])
        if rxns.annotation == []:
            rxn["annotation"] = "false"
        elif len(rxns.annotation) == 1 and  'MetaNetX (MNX) Equation' in str(rxns.annotation): 
            rxn["annotation"] = "true"
        else:
            rxn["annotation"] = "true"


def check_rxn_balance(model_info, model): 
    """
    Check if the reaction is balanced
    """ 
    for rxn in model_info['reactions']:
        if rxn['id'] not in model_info['exchange_rxns']:
            # check is a dictionary, calculate the mass and charge balance of the reaction, for a balanced reaction, this should be empty
            check = model.reactions.get_by_id(rxn['id']).check_mass_balance()
            if len(check) != 0:
                # Elemental balance, charge nan  {'charge': nan}
                if 'charge' in check.keys() and len(check) == 1 and math.isnan(check['charge']):
                    rxn["balance"] = "true"
                # {'charge': 1.0, 'H': 1.0}，The number of charge and the number of hydrogen are the same and have the same sign, the reaction is considered to be correct, and no manual inspection is performed
                elif  "'H': " in str(check) and "'charge': " in str(check) and len(check) == 2 and check['charge'] == check['H'] and check['charge'] * check['H'] >=0:
                    rxn["balance"] = "true"
                # h2o，and the same sign, this response is considered correct
                elif "'H': " in str(check) and "'O': " in str(check) and len(check) == 2 and abs(check['O']) == 1 and abs(check['H']) == 2 and check['O'] * check['H'] >=0:
                    rxn["balance"] = "true"
                else:
                    rxn["balance"] = "false"
            else:
                rxn["balance"] = "true"
        else:
            rxn["balance"] = ""


def check_C_balance(model_info, model): 
    """
    Check if the reaction is carbon balanced
    """ 
    for rxn in model_info['reactions']:
        if rxn['id'] not in model_info['exchange_rxns']:
            check = model.reactions.get_by_id(rxn['id']).check_mass_balance()
            if len(check) != 0 and 'C' in check.keys():
                # The element contains C, that is, carbon is not conserved
                rxn["c_balance"] = "false"
            else:
                rxn["c_balance"] = "true"
        else:
            rxn["c_balance"] = ""


def net_penalty_points(model_info):
    """
    give reaction penalty
    """
    # rules.get_all_rules(model_info, model)
    # check_rxn_balance(model_info, model)
    # check_C_balance(model_info, model)
    for rxn in model_info['reactions']:
        grate = 0
        if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns'] and not rxn['id'].startswith('ADD_'):
            if rxn['balance'] == 'false':  # Response Imbalance Response Penalty
                grate += 3
                rxn["rules"]["not_balance_rxn"] = "true"
            if rxn['c_balance'] == 'false':  # carbon imbalance reaction penalty
                grate += 20
                rxn["rules"]["not_c_balance_rxn"] = "true"
            if rxn["rules"]:  # 3 points for failure to comply with reaction direction rules
                grate += 3
            if rxn["annotation"] == "false":
                grate += 1
                rxn["rules"]["no_annotation_rxn"] = "true"
            rxn["net_penalty_points"] = grate

            if 'right_chain_rxn' in rxn['rules'].keys():
                rxn["net_penalty_points"] = 0
        elif rxn['id'] in model_info['exchange_rxns'] and 'o2s_rxn' in rxn["rules"].keys():
            rxn["net_penalty_points"] = 3
        else:
            rxn["net_penalty_points"] = 0
 


def write_flux_file(model_info, model, objectiveId):
    """
    get flux file
    """
    pfba_solution = pfba(model)  # Reaction with limited flux value v greater than 1e-6
    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]  # abs(pfba_solution.fluxes)>1e-6 gets the true or false result of each reaction, and finally returns true
    with open('mqc/flux.txt', 'w') as flux:  
        for r,v in need_fluxes.items():
            for rxn in model_info['reactions']:
                if r == rxn['id'] and r != objectiveId:
                    rxns = model.reactions.get_by_id(r)
                    try:
                        check = rxns.check_mass_balance()        
                        if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                            flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_name']}\t{check}\t{rxn['net_penalty_points']}\t{rxns.bounds}\n")   
                        else:
                            flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_id']}\t{check}\t{rxn['net_penalty_points']}\t{rxns.bounds}\n")   
                    except ValueError:
                        if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                            flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_name']}\t{check}\t{rxn['net_penalty_points']}\t{rxns.bounds}\n")   
                        else:
                            flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_id']}\t{rxn['net_penalty_points']}\t{rxns.bounds}\tno_check_mass_balance\n") 
                    except KeyError:
                        print(rxn, repr(KeyError))
                        exit()
            if r == objectiveId:
                rxns = model.reactions.get_by_id(r)
                if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                    flux.write(f"{r}\t{round(v,5)}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.check_mass_balance()}\t{rxns.bounds}")
                else:
                    flux.write(f"{r}\t{round(v,5)}\t{rxns.build_reaction_string()}\t{rxns.check_mass_balance()}\t{rxns.bounds}")
 

def write_flux_file2(model_info, model):
    """
    get flux file
    """
    pfba_solution = pfba(model)  # Reaction with limited flux value v greater than 1e-6
    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]  # abs(pfba_solution.fluxes)>1e-6 gets the true or false result of each reaction, and finally returns true
    with open('mqc/flux2.txt', 'w') as flux:  
        for r,v in need_fluxes.items():
            for rxn in model_info['reactions']:
                if r == rxn['id']:
                    rxns = model.reactions.get_by_id(r)
                    try:
                        check = rxns.check_mass_balance()           
                        flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_id']}\t{check}\t{rxns.bounds}\n")   
                    except ValueError:
                        flux.write(f"{r}\t{round(v,5)}\t{rxn['rxn_exp_id']}\t{rxns.bounds}\tno_check_mass_balance\n")  


def add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, identifier, num):
    """
    Add the checked reaction back to see if target product can still be generated infinitely
    """
    print('grate_delete_rxn: ',grate_delete_rxn)
    infinite_rxn = []
    for rxn in model_info['reactions']:
        # rxn[f"{identifier}_modify"] = 'false'
        for rxnId in grate_delete_rxn:
            if rxnId == rxn['id']:
                model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(rxnId).bounds
                if model.slim_optimize() >= num:
                    pfba_solution = pfba(model)
                    rxns = model.reactions.get_by_id(rxnId)            
                    if rxn['balance'] == 'false':
                        rxns.bounds = (0,0)
                        rxn['bounds'] = [0,0]
                        rxn[f"{identifier}_modify"] = "true"
                        infinite_rxn.append(rxnId)
                    else:
                        if pfba_solution.fluxes[rxnId] > 0:
                            rxns.bounds = (-1000,0)
                            rxn['bounds'] = [-1000,0]
                            rxn[f"{identifier}_modify"] = "true"
                            infinite_rxn.append(rxnId)
                        if pfba_solution.fluxes[rxnId] < 0:
                            rxns.bounds = (0,1000)
                            rxn['bounds'] = [0,1000]
                            rxn[f"{identifier}_modify"] = "true"
                            infinite_rxn.append(rxnId)
    return infinite_rxn

def get_net_infinite_info(model_info, model, model_control_info, identifier):
    """
    Get response corrections
    """
    bounds_modification = ''
    for rxn in model_info['reactions']:
        if f"{identifier}_modify" in rxn.keys() and rxn[f"{identifier}_modify"] == "true" and rxn['rules']:
            rxns = model.reactions.get_by_id(rxn['id'])
            if rxn['bounds'] == [0,0] : bounds_modification = "close"
            if rxn['bounds'] == [-1000,0] : bounds_modification = "<--"
            if rxn['bounds'] == [0,1000] : bounds_modification = "-->"
            infiniteInfo = {"反应ID" : rxn["id"],
                        "反应表达式" : rxn['rxn_exp_id'], 
                        "错误原因" : list(rxn['rules'].keys()),
                        "平衡情况" : rxns.check_mass_balance(),
                        "反应修正情况" : bounds_modification}
            if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                infiniteInfo['反应表达式'] = rxn['rxn_exp_name']
            model_control_info[f"{identifier}"]['模型修正情况'].append(infiniteInfo)


def boundary_restoration(model_info, model, identifier):
    """
    Response to changes each cycle, restoring bounds in time
    """
    for rxn in model_info['reactions']:
        if f"{identifier}_modify" in rxn.keys() and rxn[f"{identifier}_modify"] == "true":
            model.reactions.get_by_id(rxn['id']).bounds = rxn['bounds']     
            rxn[f"{identifier}_modify"] = 'false'


def get_final_net_fluxes(modelInfo, model, model_control_info):
    """"""
    final_netId_flux = {}
    for netId in model_control_info["nets"]["净物质:"]:
        with model:
            object_rxn_id = add_demand(modelInfo, model, netId)
            model.objective = object_rxn_id  
            final_netId_flux[netId] = model.slim_optimize()
    return final_netId_flux


def get_final_yield_fluxes(modelInfo, model, model_control_info):
    """"""
    final_yield_flux = {}
    for yieldId in model_control_info["yields"]["得率净物质:"]:
        with model:
            object_rxn_id = add_demand(modelInfo, model, yieldId)
            model.objective = object_rxn_id
            temp_max_rate = max_reduced_degree_rate(model, yieldId)
            final_yield_flux[yieldId] = temp_max_rate
    return final_yield_flux

def write_biomass(model_control_info, temp_list, bio_type):
    """"""
    model_control_info['biomasses'][f"{bio_type}:"].extend(temp_list)                       