# -*- coding: utf-8 -*-

from cobra.util.solver import linear_reaction_coefficients
import json 
from mqc.control.preprocessing_control import Preprocess
from mqc.utils import *


class ModelPreprocess():
    """
    Integrate the model's attributes into a dictionary.

    Attributes
    ----------
    diff_results : python.Dictionary
        The dictionary structure of memote.MemoteResult objects.
    configuration : memote.MemoteConfiguration
        A memote configuration structure.

    """
    def __init__(self) -> None:
        """
        Initialize the model dictionary.

        """
        self.model_info = {}
        
    
    def preprocess_met(self, obj):
        """
        Add preprocessing operations such as formula and charge to metabolites.

        Notes
        -----
        The name of the seed library needs to be mapped to the table, which is mapped by modelseed_id_to_name

        """
        for met in obj.model.metabolites:
            if met.id.startswith("cpd"):  # seed
                obj.add_seed_formula()
                self.modelseed_id_to_name = obj.bind_id_name()
                break 
            elif len(met.id.split("_")[-1]) == 1:  # bigg
                obj.add_bigg_formula()
                break
            elif "@" in met.id:  # meta
                obj.add_meta_formula()
                break
            elif "[" in met.id and "]" in met.id:  # virtual
                break
            elif met.id.startswith("C"):  # kegg
                obj.add_kegg_formula()
                break 
            else:
                """其他库--kegg,metacyc...bigg不好区分也可放进来,做成一个大表,这个表只包含formula就行"""


    def preprocess_rxn(self, obj):  
        """
        Additive exchange reaction without exchange reaction.

        Notes
        -----
        Find the intersection of all reaction IDs with the exchange reaction list

        """
        all_rxn_id = []
        for rxn in obj.model.reactions:
            all_rxn_id.append(rxn.id)
        if len(set(all_rxn_id) & set(H2O_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("h2o_exchange")
        if len(set(all_rxn_id) & set(NH4_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("nh4_exchange")
        if len(set(all_rxn_id) & set(PI_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("pi_exchange")
        if len(set(all_rxn_id) & set(H_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("h_exchange")
        if len(set(all_rxn_id) & set(O2_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("o2_exchange")
        if len(set(all_rxn_id) & set(CO2_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("co2_exchange")
        if len(set(all_rxn_id) & set(SO4_EXCHANGE_RXN)) == 0:
            obj.add_exchange_rxn("so4_exchange")


    def get_identifier(self, obj):
        """"""
        identifier = ""
        for met in obj.model.metabolites:
            if met.id.startswith("cpd"):  # seed
                identifier = "modelseed"
                break 
            elif len(met.id.split("_")[-1]) == 1:  # bigg
                identifier = "bigg"
                break 
            elif "@" in met.id:  # metaNetx
                identifier = "metaNetx"
                break 
            elif "[" in met.id and "]" in met.id:  # virtual
                identifier = "virtual"
                break
            else:
                identifier = "other"
                break 
        return identifier


    def get_met_info(self, obj):
        """
        Collect metabolite information.

        Returns
        -------
        met_info : list
            metabolite information

        """
        met_info = []
        for met in obj.model.metabolites:
            if met.id.startswith("cpd"):  # seed
                metInfo = {"id" : met.id,
                            "name" : self.modelseed_id_to_name[met.id.split('_')[0]],
                            "charge" : met.charge,
                            "formula" : met.formula}
            elif len(met.id.split("_")[-1]) == 1:  # bigg
                metInfo = {"id" : met.id,
                            "name" : '_'.join(met.id.split('_')[:-1]),
                            "charge" : met.charge,
                            "formula" : met.formula} 
            elif "@" in met.id:  # metaNetx
                metInfo = {"id" : met.id,
                            "name" : met.name,
                            "charge" : met.charge,
                            "formula" : met.formula}
            elif "[" in met.id and "]" in met.id:  # virtual
                metInfo = {"id" : met.id,
                            "name" : met.id.split('[')[0],
                            "charge" : met.charge,
                            "formula" : met.formula}
            else:  # IDs that do not belong to the four libraries
                metInfo = {"id" : met.id,  
                            "name" : met.name,
                            "charge" : met.charge,
                            "formula" : met.formula}
            met_info.append(metInfo)
        return met_info

    def get_met_name(self, obj, metId):
        met = obj.model.metabolites.get_by_id(metId)
        if metId.startswith("cpd") : return self.modelseed_id_to_name[metId.split('_')[0]]  # seed
        elif len(metId.split("_")[-1]) == 1 : return '_'.join(metId.split('_')[:-1])  # bigg 
        elif "@" in metId : return met.name  # metaNetx
        elif "[" in metId and "]" in met.id : return metId.split('[')[0]  # virtual
        else : return met.name  # IDs that do not belong to the four libraries

    def get_rxn_info(self, obj):
        """
        Collect reaction information.

        Returns
        -------
        rxn_info : list
            reaction information

        """
        rxn_info = []
        for rxn in obj.model.reactions:
            rxnInfo = {"id" : rxn.id,
                        "name" : rxn.name,
                        "bounds" : rxn.bounds,
                        # "lower_bound" : rxn.lower_bound,
                        # "upper_bound" : rxn.upper_bound,
                        "rxn_exp_id" : rxn.build_reaction_string(),
                        "rxn_exp_name" : rxn.build_reaction_string(use_metabolite_names = True),
                        "annotation" : rxn.annotation}
            rxnInfo["all_mets"] = [self.get_met_name(obj, met.id) for met in rxn.metabolites]   
            rxnInfo["reactants_mets"] = [self.get_met_name(obj, met.id) for met in rxn.reactants] 
            rxnInfo["products_mets"] = [self.get_met_name(obj, met.id) for met in rxn.products] 
            rxnInfo["rules"] = {}
            rxn_info.append(rxnInfo)
        return rxn_info


    def get_exchange_rxns(self):
        """
        Get exchange reaction ID.

        Returns
        -------
        exchange_rxns : list
            exchange reaction id

        """
        exchange_rxns = []
        for rxn in self.model_info["reactions"]:
            if len(rxn["all_mets"]) == 1:
                exchange_rxns.append(rxn["id"])
        return exchange_rxns

    def get_transport_rxns(self):
        """
        Get the transport reaction ID.

        Returns
        -------
        transport_rxns : list
            transport reaction id

        """
        transport_rxns = []
        for rxn in self.model_info["reactions"]:
            if sorted(rxn["reactants_mets"]) == sorted(rxn["products_mets"]):
                transport_rxns.append(rxn["id"])
        return transport_rxns


    def get_initial_obj_info(self, obj):

        # initial_rxn_info = []
        try:
            initial_rxn_id = list(linear_reaction_coefficients(obj.model).keys())[0].id
        except:
            initial_rxn_id = set_initial_obj(obj.model)
        else:
             initial_rxn_id = set_initial_obj2(obj.model, initial_rxn_id)
        initial_rxn_flux = obj.model.slim_optimize()
        try:
            if self.model_info['model_identifier'] == 'modelseed' or self.model_info['model_identifier'] == 'metaNetx':
                initial_rxn_exp = obj.model.reactions.get_by_id(initial_rxn_id).build_reaction_string(use_metabolite_names=True)
            else:
                initial_rxn_exp = obj.model.reactions.get_by_id(initial_rxn_id).build_reaction_string()
        except KeyError:
            initial_rxn_exp = ''
        InitialRxnInfo = {"initial_rxn_id" : initial_rxn_id,
                            "initial_rxn_flux" : initial_rxn_flux,
                            "initial_rxn_exp" : initial_rxn_exp}
        # initial_rxn_info.append(InitialRxnInfo)
        return InitialRxnInfo

    def find_glucose(self):
        """"""
        glucose_rxnId = ''
        if len(set(GLUCOSE_RXN) & set(self.model_info['all_rxn_id'])) != 0:
            glucose_rxnId = list(set(GLUCOSE_RXN) & set(self.model_info['all_rxn_id']))[0]
        else:
            for rxn in self.model_info["reactions"]:
                if rxn['id'] in self.model_info['exchange_rxns'] and 'glucose' in rxn['all_mets'][0].lower():
                    glucose_rxnId = rxn['id']
        return glucose_rxnId


    def get_model_info(self, controler):
        """
        Integrate model properties together.

        Parameters
        ----------
        file_path : str
            model file

        """
        self.preprocess_met(controler)
        self.preprocess_rxn(controler)
  
        self.model_info = {"model_id" : controler.model.id,
                           "model_identifier" : self.get_identifier(controler),
                            "metabolites" : self.get_met_info(controler),
                            "reactions" : self.get_rxn_info(controler),}
        self.model_info["initial_rxn"] = self.get_initial_obj_info(controler)
        self.model_info["all_rxn_id"] = [rxn.id for rxn in controler.model.reactions]
        self.model_info["exchange_rxns"] = self.get_exchange_rxns()
        self.model_info["transport_rxns"] = self.get_transport_rxns()
        self.model_info["glucose_rxnId"] = self.find_glucose()
   
        modelInfo = json.dumps(self.model_info, ensure_ascii=False, allow_nan = True, indent=1)
        with open("mqc/test.json", "w", newline='',) as f:
            f.write(modelInfo)



