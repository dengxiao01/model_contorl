# -*- coding: utf-8 -*-

import re 
from mqc.defaults import *

class Rules():

    def __init__(self):


    def find_o2s(model_info):
        """
        o2s rules
        """
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:    
                if rxn["all_mets"][0] in O2S:
                    rxn["rules"]["o2s_rxn"] = "true"

    def find_O2_inorganic_rxn(model_info):
        """
        find inorganic reactions involving oxygen
        """
        O2_inorganic_rxn = []   
        for rxn in model_info["reactions"]:
            if rxn["id"] not in model_info["transport_rxns"] and rxn["id"] not in model_info["exchange_rxns"]:
                n = 0
                s = 0
                for met in rxn["all_mets"]:
                    s += 1
                    if met in O2_METS:
                        n += 1
                if s == n:
                    O2_inorganic_rxn.append(rxn["id"])
        return O2_inorganic_rxn


    def find_superoxide_rxn(model_info):
        """
        superoxide, reaction of peroxide to produce oxygen

        Notes
        -----
        SPODMpp: 2.0 h_p + 2.0 o2s_p --> h2o2_p + o2_p
        SOD_m: 2.0 sox_m --> h2o2_m + o2_m
        GTHP_CAT: 2.0 gthrd_c + 3.0 h2o2_c --> gthox_c + 4.0 h2o_c + o2_c   Reaction of peroxide and reduced glutathione to produce oxygen

        """
        superoxide_rxn = []
        for rxn in model_info["reactions"]:
            if len(set(O2S) & set(rxn['reactants_mets'])) != 0 or len(set(SOX) & set(rxn['reactants_mets'])) != 0:
                if len(set(H2O2) & set(rxn['products_mets'])) != 0 and len(set(O2_NAME) & set(rxn['products_mets'])) != 0:
                    superoxide_rxn.append(rxn['id'])
            if len(set(O2S) & set(rxn['products_mets'])) != 0 or len(set(SOX) & set(rxn['products_mets'])) != 0:
                if len(set(H2O2) & set(rxn['reactants_mets'])) != 0 and len(set(O2_NAME) & set(rxn['reactants_mets'])) != 0:
                    superoxide_rxn.append(rxn['id'])
            if len(set(GTHRD) & set(rxn['reactants_mets'])) != 0 or len(set(H2O2) & set(rxn['reactants_mets'])) != 0:
                if len(set(GTHOX) & set(rxn['products_mets'])) != 0 and len(set(O2_NAME) & set(rxn['products_mets'])) != 0:
                    superoxide_rxn.append(rxn['id'])       
            if len(set(GTHRD) & set(rxn['products_mets'])) != 0 or len(set(H2O2) & set(rxn['products_mets'])) != 0:
                if len(set(GTHOX) & set(rxn['reactants_mets'])) != 0 and len(set(O2_NAME) & set(rxn['reactants_mets'])) != 0:
                    superoxide_rxn.append(rxn['id'])
        return superoxide_rxn


    def find_photosynthetic_rxn(model_info):
        """
        What are the photosynthetic bacteria, the reaction of photosynthesis to produce oxygen
        """
        photosynthetic_rxn = []
        for rxn in model_info['reactions']:
            if re.search(r'Photosystem',rxn['name']):
                photosynthetic_rxn.append(rxn['id'])
        return photosynthetic_rxn

    def find_O2_rxn(self, model_info):
        """
        Find the reaction that obeys the oxygen rule
        """
        O2_rxn = []
        O2_inorganic_rxn = self.find_O2_inorganic_rxn(model_info)
        superoxide_rxn = self.find_superoxide_rxn(model_info)
        photosynthetic_rxn = self.find_photosynthetic_rxn(model_info)
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns'] and rxn['id'] not in O2_inorganic_rxn and rxn['id'] not in superoxide_rxn and rxn['id'] not in photosynthetic_rxn:
                if len(set(O2_NAME) & set(rxn['reactants_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [-1000,0]:
                        rxn["rules"]["o2_rxn"] = "true"
                        O2_rxn.append(rxn['id'])
                if len(set(O2_NAME) & set(rxn['products_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [0,1000]:
                        rxn["rules"]["o2_rxn"] = "true"
                        O2_rxn.append(rxn['id'])
        return O2_rxn
        
    def find_nh4_inorganic_rxn(model_info):
        """
        find the inorganic reaction of nh3 and nh4
        """
        nh4_inorganic_rxn=[] 
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns']:
                if len(set(NH3_NAME) & set(rxn['reactants_mets'])) != 0 and len(set(NH4_NAME) & set(rxn['products_mets'])) != 0:
                    nh4_inorganic_rxn.append(rxn['id'])
                if len(set(NH3_NAME) & set(rxn['products_mets'])) != 0 and len(set(NH4_NAME) & set(rxn['reactants_mets'])) != 0:
                    nh4_inorganic_rxn.append(rxn['id'])
        return nh4_inorganic_rxn

    def find_nh3_nh4_rxn(self, model_info):
        """
        Except nh3, nh4 and ATP, nadh, nadph, chor, prpp, udpacblfuc reactions can fix ammonia, others are reactions to generate ammonia
        """
        nh4_inorganic_rxn = self.find_nh4_inorganic_rxn(model_info)
        O2_rxn = self.find_O2_rxn(model_info)
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns'] and rxn['id'] not in nh4_inorganic_rxn and rxn['id'] not in O2_rxn:
                if len(set(NH4_NAME) & set(rxn['reactants_mets'])) != 0 or len(set(NH3_NAME) & set(rxn['reactants_mets'])) != 0:
                    if len(set(SOLID_AMMONIA) & set(rxn['reactants_mets'])) == 0:
                        if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [0,1000]:                            
                            rxn["rules"]["nh3_nh4_rxn"] = "true"                 
                if len(set(NH4_NAME) & set(rxn['products_mets'])) != 0 or len(set(NH3_NAME) & set(rxn['products_mets'])) != 0:
                    if  len(set(SOLID_AMMONIA) & set(rxn['products_mets'])) == 0:
                        if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [-1000,0]:
                            rxn["rules"]["nh3_nh4_rxn"] = "true"


    def find_CO2_rxn_all(model_info):
        """
        Find all CO2 reactions
        """
        CO2_rxn_all = []
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns']:
                if rxn['bounds'] == [-1000,1000] and len(set(CO2_NAME) & set(rxn['all_mets'])) != 0:
                    CO2_rxn_all.append(rxn['id'])
                if rxn['bounds'] == [0,1000] and len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0:
                    CO2_rxn_all.append(rxn.id)
                if rxn['bounds'] == [-1000,0] and len(set(CO2_NAME) & set(rxn['products_mets'])) != 0:
                    CO2_rxn_all.append(rxn.id)
        return CO2_rxn_all     

    def find_co2_inorganic_rxn(self, model_info):
        """
        find the inorganic reaction involving co2
        """
        co2_inorganic_rxn=[]
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        for rxn in model_info['reactions']:
            if rxn['id'] in CO2_rxn_all and set(rxn['all_mets']).issubset(CO2_INORGANIC):
                co2_inorganic_rxn.append(rxn['id'])
        return co2_inorganic_rxn   

    def find_co2_atp_pep_rxn(self, model_info):
        """
        co2 and high energy compound atp pep can fix carbon
        """
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_atp_pep_rxn=[]
        for rxn in model_info['reactions']:
            if rxn['id'] in CO2_rxn_all:
                if len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0:
                    if len(set(ATP_NAME) & set(rxn['reactants_mets'])) != 0 or len(set(PEP) & set(rxn['reactants_mets'])) != 0:
                        co2_atp_pep_rxn.append(rxn['id'])
                
                if len(set(CO2_NAME) & set(rxn['products_mets'])) != 0:
                    if len(set(ATP_NAME) & set(rxn['products_mets'])) != 0 or len(set(PEP) & set(rxn['products_mets'])) != 0:
                        co2_atp_pep_rxn.append(rxn['id'])
        return co2_atp_pep_rxn   

    def find_natural_co2_rxn(self, model_info):
        """
        6 natural carbon fixation reactions
        """
        natural_co2_rxn = []
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_inorganic_rxn = self.find_co2_inorganic_rxn(model_info)
        # Calvin cycle RBPC: co2_c + h2o_c + rb15bp_c --> 2.0 3pg_c + 2.0 h_c
        for rxn in model_info['reactions']:
            if rxn['id'] in CO2_rxn_all and rxn['id'] not in co2_inorganic_rxn:
                if len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0 and len(set(RB15BP) & set(rxn['reactants_mets'])) != 0 and len(set(PG3) & set(rxn['products_mets'])) != 0:
                    natural_co2_rxn.append(rxn['id'])
                if len(set(CO2_NAME) & set(rxn['products_mets'])) != 0 and len(set(RB15BP) & set(rxn['products_mets'])) != 0 and len(set(PG3) & set(rxn['reactants_mets'])) != 0:
                    natural_co2_rxn.append(rxn['id'])
                # for carbon fixation reaction
                if len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0 and len(set(FOR) & set(rxn['products_mets'])) != 0: 
                    natural_co2_rxn.append(rxn['id'])
                if len(set(CO2_NAME) & set(rxn['products_mets'])) != 0 and len(set(FOR) & set(rxn['reactants_mets'])) != 0:
                    natural_co2_rxn.append(rxn['id'])  
                
                if len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0 and len(set(FOR) & set(rxn['products_mets'])) != 0: 
                    natural_co2_rxn.append(rxn['id'])
                if len(set(CO2_NAME) & set(rxn['products_mets'])) != 0 and len(set(FOR) & set(rxn['reactants_mets'])) != 0:
                    natural_co2_rxn.append(rxn['id']) 
                # CODH_ACS: co2_c + coa_c + fdxr_42_c + h_c + mecfsp_c --> accoa_c + cfesp_c + fdxo_42_c + h2o_c
                natural_co2_rxn.append('CODH_ACS')
                
                # reverse TCA cycle
                # OOR2r: akg_c + coa_c + fdxo_42_c <=> co2_c + fdxr_42_c + h_c + succoa_c
                # ICDHyr: icit_c + nadp_c <=> akg_c + co2_c + nadph_c
                if set(CO2_TCA_I).issubset(rxn['all_mets']): 
                    natural_co2_rxn.append(rxn['id'])
                if set(CO2_TCA_II).issubset(rxn['all_mets']): 
                    natural_co2_rxn.append(rxn['id'])
                # ferredoxin, pyr to generate accoa
                # POR5: coa_c + 2.0 flxso_c + pyr_c <=> accoa_c + co2_c + 2.0 flxr_c + h_c
                # POR: coa_c + fdxo_42_c + pyr_c <=> accoa_c + co2_c + fdxr_42_c + h_c
                if set(CO2_ACCOA_I).issubset(rxn['all_mets']): 
                    natural_co2_rxn.append(rxn['id'])
                if set(CO2_ACCOA_II).issubset(rxn['all_mets']): 
                    natural_co2_rxn.append(rxn['id'])
        return natural_co2_rxn

    def find_CO2_rxn(self, model_info):
        """
        Find reactions that do not comply with the CO2 rule
        """
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_inorganic_rxn = self.find_co2_inorganic_rxn(model_info)
        co2_atp_pep_rxn = self.find_co2_atp_pep_rxn(model_info)
        natural_co2_rxn = self.find_natural_co2_rxn(model_info)
        for rxn in model_info['reactions']:
            if rxn['id'] in list(set(CO2_rxn_all) - set(co2_inorganic_rxn) - set(co2_atp_pep_rxn) -set(natural_co2_rxn)):
                if len(set(CO2_NAME) & set(rxn['reactants_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [0,1000]:
                        rxn["rules"]["co2_rxn"] = "true"
                if len(set(CO2_NAME) & set(rxn['products_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [-1000,0]:
                        rxn["rules"]["co2_rxn"] = "true"


    def find_ATP_synthase_rxn(model_info):
        """
        Find out the list of ATP synthase reactions
        """
        ATP_synthase_rxn=[]
        for rxn in model_info['reactions']: 
            if sorted(rxn['all_mets']) == sorted(ATP_SYNTHASE):
                ATP_synthase_rxn.append(rxn['id'])
        return ATP_synthase_rxn                

    def find_ATP_rxn(self, model_info):
        """
        Find reactions that do not comply with ATP rules
        """
        ATP_synthase_rxn =  self.find_ATP_synthase_rxn(model_info)
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns'] and rxn['id'] not in ATP_synthase_rxn:
                if len(set(XTP) & set(rxn['reactants_mets'])) != 0:
                    if len(set(XPI) & set(rxn['products_mets'])) != 0:
                        if len(set(rxn['products_mets']) & set(YLCOA)) == 0 and rxn['bounds'][0] == -1000:  
                            rxn["rules"]["atp_rxn"] = "true"
                if len(set(XTP) & set(rxn['products_mets'])) != 0:
                    if len(set(XPI) & set(rxn['reactants_mets'])) != 0:
                        if len(set(rxn['reactants_mets']) & set(YLCOA)) == 0 and rxn['bounds'][0] == -1000:   
                            rxn["rules"]["atp_rxn"] = "true"


    def find_proton_rxn(self, model_info):
        """
        find the proton transport reaction
        """
        proton_rxn, atp_synthetic = [], []
        hm = []
        for rxn in model_info['reactions']:

            h_met = [met.id for met in rxn['all_mets'] if met == 'h']
            if sorted(ATP_SYNTHASE) == sorted(rxn_met):
                if 'atp' in r_mets :
                    if r.bounds == (0,1000) : proton_rxn.append(r.id)
                    h_met.reverse()
                    hm = h_met
                if 'atp' in pro_mets :
                    if r.bounds == (-1000,0) : proton_rxn.append(r.id)
                    hm = h_met
            
                for rxn in model.reactions:
                    rxn_met_c = [met.id for met in rxn.metabolites]
                    r_mets_id = [str(m) for m in rxn.reactants]
                    pro_mets_id = [str(m) for m in rxn.products]
                    if len(set(hm)&set(rxn_met_c)) == 2 and rxn.id in transport_rxns:
                        if hm[0] in r_mets_id and hm[1] in pro_mets_id and rxn.bounds == (-1000,0):
                            proton_rxn.append(rxn.id)
                        elif hm[0] in pro_mets_id and hm[1] in r_mets_id and rxn.bounds == (0,1000):            
                            proton_rxn.append(rxn.id)
                hm = []
        return proton_rxn, atp_synthetic                        

    def find_Sugar_hydrolysis_rxn(model_info):
        """
        Find reactions that don't follow the rules for polysaccharides and glycogen
        """
        sugar = []
        for sugar_met in SUGAR:
            for met in model_info['metabolites']:
                if sugar_met in met['name']:
                    sugar.append(met['name'])

        sugar = list(set(sugar)- set(EXCEPT_SUGAR))

        # polysaccharide hydrolysis
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns'] and rxn['id'] not in model_info['transport_rxns']:
                if len(set(sugar) & set(rxn['reactants_mets'])) == 1 and len(set(H2O) & set(rxn['reactants_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [-1000,0]:
                        rxn["rules"]["sugar_hydrolysis_rxn"] = "true"
                if len(set(sugar) & set(rxn['products_mets'])) == 1 and len(set(H2O) & set(rxn['products_mets'])) != 0:
                    if rxn['bounds'] == [-1000,1000] or rxn['bounds'] == [0,1000]:
                        rxn["rules"]["sugar_hydrolysis_rxn"] = "true"

    def find_PTS_transport(model_info):
        """
        Find reactions that do not correspond to the PTS pathway for material transport. PTS pathway for material transport, irreversible
        """
        for rxn in model_info['reactions']:
            if 'PTS' in rxn['name'] and rxn['bounds'] == [-1000,1000]:
                rxn["rules"]["PTS_transport"] = "true"
                
    def find_ppi_h2o(model_info):
        """
        Hydrolysis of polyphosphoric acid to generate pi ppi
        """
        for rxn in model_info['reactions']:
            if rxn['id'] not in model_info['exchange_rxns']:
                if len(set(H2O) & set(rxn['reactants_mets'])) != 0:
                    for r_m in rxn['reactants_mets']:
                        if re.search(r'triphosphate|diphosphate',r_m) and rxn['bounds'][0] == -1000: 
                            if len(set(PI_NAME) & set(rxn['products_mets'])) != 0 or len(set(PPI_NAME) & set(rxn['products_mets'])) != 0:
                                rxn["rules"]["ppi_h2o"] = "true"       
                if len(set(H2O) & set(rxn['products_mets'])) != 0:            
                    for p_m in rxn['products_mets']:
                        if re.search(r'triphosphate|diphosphate',p_m) and rxn['bounds'][1] == 1000: 
                            if len(set(PI_NAME) & set(rxn['reactants_mets'])) != 0 or len(set(PPI_NAME) & set(rxn['reactants_mets'])) != 0:
                                rxn["rules"]["ppi_h2o"] = "true"             