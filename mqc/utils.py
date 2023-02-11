# -*- coding: utf-8 -*-

from cobra import Model, Reaction, Metabolite
import pandas as pd
import math 
from mqc.defaults import *

def get_model_file():
    """
    Get the model file.

    Return
    ------
    str
        model file

    """
    file = "mqc/test_data/virtual/Acinetobacter_haemolyticus_NIPH_261.xml"
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
        all_mets = rxn["all_mets"]
        if len(all_mets) == 1 and metId in all_mets:
            objectiveId = rxn["id"]
            flag += 1
    if flag == 0:
        if 'DM_' + metId in modelInfo["all_rxn_obj"]:
            objectiveId = f"DM_metId"
        else:
            demand_rxn = model.add_boundary(mets, type = 'demand') # 添加demand反应
            objectiveId = demand_rxn.id
    return objectiveId


def set_model_objective(modelInfo, model, rxnId):
    """
    Set the incoming response ID as the target and update modelInfo.
    """
    now_obj_rxn_info = []
    model.objective = rxnId
    now_obj_rxn_flux = model.optimize().objective_value
    now_obj_rxn_exp = model.reactions.get_by_id(rxnId).build_reaction_string()   
    NowObjRxnInfo = {"now_obj_rxn_id" : rxnId,
                    "now_obj_rxn_flux" : now_obj_rxn_flux,
                    "now_obj_rxn_exp" : now_obj_rxn_exp}
    now_obj_rxn_info.append(NowObjRxnInfo)
    modelInfo["now_obj_rxn"] = now_obj_rxn_info


def is_bio_exp_correct(modelInfo):
    """
    Judging whether the biomass response of the model is correctly represented.
    """
    id_index = modelInfo["initial_rxn"].index("initial_rxn_id")
    if not modelInfo["initial_rxn"][id_index]:
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

def check_rxn_balance(model_info, model, rxnId): 
    """
    Check if the reaction is balanced
    """ 
    for rxn in model_info['reactions']:
        if rxn['id'] == rxnId and rxnId not in model_info['exchange_rxns']:
            # check is a dictionary, calculate the mass and charge balance of the reaction, for a balanced reaction, this should be empty
            check = model.reactions.get_by_id(rxnId).check_mass_balance()
            if len(check) != 0:
                # Elemental balance, charge nan  {'charge': nan}
                if 'charge' in check.keys() and len(check) == 1 and math.isnan(check['charge']):
                    rxn["balance"] = "true"
                # {'charge': 1.0, 'H': 1.0}，The number of charge and the number of hydrogen are the same and have the same sign, the reaction is considered to be correct, and no manual inspection is performed
                elif  "'H': " in str(check) and "'charge': " in str(check) and len(check) == 2 and check['charge'] == check['H'] and check['charge'] * check['H'] >=0:
                    rxn["balance"] = "true"
                # h2o，and the same sign, this response is considered correct
                elif "'H': " in str(check) and "'O': " in str(check) and len(check) == 2 and abs(check['O']) == 1 and  abs(check['H']) == 2 and check['O'] * check['H'] >=0:
                    rxn["balance"] = "true"
                else:
                    rxn["balance"] = "false"
            rxn["balance"] = "true"

def check_C_balance(model_info, model, rxnId): 
    """
    Check if the reaction is carbon balanced
    """ 
    for rxn in model_info['reactions']:
        if rxn['id'] == rxnId and rxnId not in model_info['exchange_rxns']:
            check = model.reactions.get_by_id(rxnId).check_mass_balance()
            if len(check) != 0 and 'C' in check.keys():
                # The element contains C, that is, carbon is not conserved
                rxn["c_balance"] = "false"
            rxn["c_balance"] = "true"