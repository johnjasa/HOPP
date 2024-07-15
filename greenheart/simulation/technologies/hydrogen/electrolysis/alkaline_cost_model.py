import pandas as pd
import numpy as np

class AlkalineCosts:
    def __init__(self,alkaline_design,plant_life_yrs:int,operation_start_year:int):
        self.plant_life = plant_life_yrs
        self.operation_start_year = operation_start_year
        years_of_operation = np.arange(int(self.operation_start_year),int(self.operation_start_year+self.plant_life),1)
        self.years_of_operation = [str(y) for y in years_of_operation]
        self.pf_dict = {}
        pass
    def baseline_alkaline_costs(self):
        stack_replacement_cost_percentage = 0.15
        #cheaper capex in $/kW than PEM
        uninstalled_capex = [] #less than PEM
        #More maintenance compared to PEM
        variable_om = [] #higher than PEM
        fixed_om = [] #higher than PEM
        #Built on-site and larger footprint than PEM
        installation_factor = [] #higher than PEM
        indirect_percentage = []
        water_feedstock_cost = [] #$/gal-H2O
        KOH_feedstock_cost = [] #$/gram KOH

        pass
    def user_specified_costs(self):
        pass
    def add_profast_params(self):
        
        pf_params_dict = {}
        pf_params_dict.update({"commodity":{"initial price":100,"name":"Hydrogen","unit":"kg","escalation":0.0}})
        pf_params_dict.update({"capacity":nameplate_kg_pr_day})
        pf_params_dict.update({"long term utilization":capacity_factor})
        pf_params_dict.update({"operating life":self.plant_life})
        self.pf_dict.update({"params":pf_params_dict})

    def add_profast_capital_items(self):
        self.pf_dict.update({"capital_items":pf_capex_dict})
        pass
    def add_profast_fixed_costs(self):
        self.pf_dict.update({"fixed_costs":pf_fopex_dict})
        pass
    def add_profast_feedstock_costs(self):
        self.pf_dict.update({"feedstocks":pf_fopex_dict})
        pass
    


    def create_alkaline_profast_info(self,H2_Results):
        if "Performance By Year" in H2_Results.keys():
            lta = H2_Results["Performance By Year"]
            lta.index = self.years_of_operation
        