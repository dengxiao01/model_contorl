import cobra
from cobra.flux_analysis import pfba
import pandas as pd
import math
import re
import os
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import xlrd3
import re
import openpyxl
from multiprocessing.pool import Pool
from multiprocessing import Process, Lock
import time 


def add_rxn(rxn_name, rxn_expression, model):
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


def add_bigg_formula(model):
    """
    Add formula and charge to bigg.

    Parameters
    ----------
    model : cobra.Model

    Return
    ------
    cobra.Model
        Model after adding formula and charge

    """
    bigg_pd = pd.read_excel('workspace/model_quality_control/bigg-met.xlsx')
    bigg_dict_formula = dict(zip(bigg_pd['id'], bigg_pd['formula']))
    bigg_dict_charge = dict(zip(bigg_pd['id'], bigg_pd['charge']))
    for met in model.metabolites:
        if met.id in list(bigg_pd['id']):
            met.formula = bigg_dict_formula[f"{met.id}"]
            met.charge = bigg_dict_charge[f"{met.id}"]
    return model


def add_meta_formula(model):  
    """
    Add formula to meta.

    Parameters
    ----------
    model : cobra.Model

    Return
    ------
    cobra.Model
        Model after adding formula 

    """
    meta_met = pd.read_excel('/home/dengxiao/workspace/cobrapy/unify_id_name/meta_id_names7.xlsx')
    meta_formula = dict(zip(meta_met['meta_id'], meta_met['formula']))
    for met in model.metabolites:
        if met.formula:
            pass
        else:
            if str(met.id).split('@')[0] in list(meta_formula.keys()):
                met.formula = meta_formula[str(met.id).split('@')[0]]
    return model
    


def add_seed_formula(model):
    """
    Add formula to modelseed.

    Parameters
    ----------
    model : cobra.Model

    Returns
    -------
    cobra.Model
        Model after adding formula 

    """
    df = pd.read_excel('workspace/cobrapy/modelseed_del/modelseed_metabolites.xlsx')
    modelseed_id_to_name = dict(zip(df['id'], df['name']))
    modelseed_formula = dict(zip(df['id'], df['formula']))
    for met in model.metabolites:
        if met.id.startswith('cpd'):
            if met.formula:
                pass
            else:
                met.formula = str(modelseed_formula[met.id.split('_')[0]])
    return model, modelseed_id_to_name


def bind_id_name():
    """
    Through the original form, bind the reaction id and name.

    Parameters
    ----------
    model : cobra.Model

    Returns
    -------
    cobra.Model
        Model after adding formula 
    
    Notes
    -----
    The name directly read from the model is very messy, and it is difficult to directly extract the name

    """
    df = pd.read_excel('workspace/cobrapy/modelseed_del/modelseed_metabolites.xlsx')
    modelseed_id_to_name = dict(zip(df['id'], df['name']))
    return modelseed_id_to_name