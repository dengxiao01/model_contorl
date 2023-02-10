# -*- coding: utf-8 -*-

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
from mqc.defaults import *
from mqc.utils import add_rxn

class Preprocess():
    """
    Model preprocessing related.

    Attributes
    ----------
    diff_results : python.Dictionary
        The dictionary structure of memote.MemoteResult objects.
    configuration : memote.MemoteConfiguration
        A memote configuration structure.

    """
    
    def __init__(self, file_path):
        """
        load initial model.

        Parameters
        ----------
        file_path : str
            model file

        """
        try:
            self.model = cobra.io.read_sbml_model(file_path)
            self.check_model = cobra.io.read_sbml_model(file_path)
        except:
            self.model = cobra.io.load_json_model(file_path)
            self.check_model = cobra.io.load_json_model(file_path)



    def add_bigg_formula(self):
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
        bigg_pd = pd.read_excel("mqc/summary/bigg-met.xlsx")
        bigg_pd_formula = dict(zip(bigg_pd['id'], bigg_pd['formula']))
        bigg_pd_charge = dict(zip(bigg_pd['id'], bigg_pd['charge']))
        for met in self.model.metabolites:
            if met.id in list(bigg_pd['id']):
                met.formula = bigg_pd_formula[f"{met.id}"]
                met.charge = bigg_pd_charge[f"{met.id}"]


    def add_meta_formula(self):  
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
        meta_met = pd.read_excel("mqc/summary/meta_id_names.xlsx")
        meta_formula = dict(zip(meta_met['meta_id'], meta_met['formula']))
        for met in self.model.metabolites:
            if met.formula:
                pass
            else:
                if met.id.split('@')[0] in list(meta_formula.keys()):
                    met.formula = meta_formula[met.id.split('@')[0]]
        

    def add_exchange_rxn(self, met_exchange):
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
        h_c, h2o_c, pi_c, nh4_c, o2_c, co2_c, so4_c = '', '', '', '', '', '', ''
        for met in self.model.metabolites:
            if met.id in H : h_c = met.id 
            if met.id in H2O : h2o_c = met.id 
            if met.id in PI : pi_c = met.id 
            if met.id in NH4 : nh4_c = met.id 
            if met.id in O2 : o2_c = met.id 
            if met.id in CO2 : co2_c = met.id 
            if met.id in SO4 : so4_c = met.id 
        rxn_name = ['h2o_exchange', 'nh4_exchange', 'pi_exchange', 'h_exchange', 'o2_exchange', 'co2_exchange', 'so4_exchange']
        rxn_exp = [h2o_c + ' <=>', nh4_c + ' <=>', pi_c + ' <=>', h_c + ' <=>', o2_c + ' <=>', co2_c + ' <=>', so4_c + ' <=>']
        for i in range(7):
            if rxn_name[i] == met_exchange:
                add_rxn(self.model, rxn_name[i], rxn_exp[i])
                add_rxn(self.check_model, rxn_name[i], rxn_exp[i])


    def add_seed_formula(self):
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
        df = pd.read_excel("mqc/summary/modelseed_metabolites.xlsx")
        modelseed_formula = dict(zip(df['id'], df['formula']))
        for met in self.model.metabolites:
            if met.formula:
                pass
            else:
                met.formula = str(modelseed_formula[met.id.split('_')[0]])


    def bind_id_name(self):
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
        df = pd.read_excel("mqc/summary/modelseed_metabolites.xlsx")
        modelseed_id_to_name = dict(zip(df['id'], df['name']))
        return modelseed_id_to_name

