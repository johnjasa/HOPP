from greenheart.simulation.technologies.hydrogen.electrolysis.electrolyzer_base_class import ElectrolyzerCluster
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



def estimate_alkaline_footprint_singlitico(electrolyzer_size_MW):
    #https://doi.org/10.1016/j.rset.2021.100005
    footprint_m2_pr_MW = 95
    footprint_m2 = footprint_m2_pr_MW*electrolyzer_size_MW
    return footprint_m2

def estimate_alkaline_system_footprint(n_stacks):
    #Thomas I. Valdez, ... Sebastian Freund, in Machinery and Energy Systems for the Hydrogen Economy, 2022
    #https://www.sciencedirect.com/topics/engineering/alkaline-electrolysis
    # 5.4.7.2: Commercial alkaline electrolysis systems
    # 0.5 MW electrolyzer: 2.2m x 2.3 m by 3m tall
    #Teledyne Titan EL water electrolysis system in 80 Nm3/h configuration.
    #assuming double the volume of 0.5 MW electrolyzer and square shaped
    ref_stack_width = 5.06 #m
    ref_stack_length = 5.06 #m
    ref_stack_height = 3 #m
    #https://www.sinohyenergy.com/5mw-10mw-alkaline-water-electrolysis/
    #1800 L H2O/hour, 2000 Nm^3-H2/hr, 1000 Nm^3 O2/hr
    system_width = n_stacks*ref_stack_width
    system_length = n_stacks*ref_stack_length
    system_area = system_width*system_length
    system_volume = system_area*ref_stack_height
    return {"System Width [m]":system_width,"System Length [m]":system_length,"System Area [m^2]":system_area,"System Volume [m^3]":system_volume}

def gibbs_temp_pressure(T_c):
    #https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2
    #^ shomate equatiom
    T_k=T_c+ 273.15  # convert Celsius to Kelvin
    t = T_k/1000
    enthalpy_stc = -285.83 #kJ/mol
    entropy = 69.95 #J/mol-K
    a = -203.6060
    b = 1523.290
    c = -3196.413
    d = 2474.455
    e = 3.855326
    f = -256.5478
    g = -488.7163
    h = -285.8304

    delta_h = a*t + b*(t**2)/2 + c*(t**3)/3 + d*(t**4)/4 - (e/t) + f - h
    s = a*np.log(t) + b*t + c*(t**2)/2 + d*(t**3)/3 - e/(2*(t**2)) + g

    pass
def get_efficiency_curve(electrolyzer: ElectrolyzerCluster, file_desc = "test"):
    filepath = os.path.join(os.path.dirname(__file__), f"{electrolyzer.electrolyzer_type}_Efficiency-Curve-{file_desc}.csv")
    dA = 10
    current_range = np.arange(dA,electrolyzer.nominal_current+dA,dA) 
    current_density = np.zeros(len(current_range))
    V_cell = np.zeros(len(current_range))
    H2_cell = np.zeros(len(current_range)) #kg/hr

    for ii,I_stack in enumerate(current_range):
        if electrolyzer.electrolyzer_type == "ALK":
            current_density[ii] = electrolyzer.calc_current_density(I_stack, electrolyzer.T_stack)
        else:
            current_density[ii] = electrolyzer.calc_current_density(I_stack)

        V_cell[ii] = electrolyzer.cell_design(electrolyzer.T_stack,I_stack)
        
        if electrolyzer.electrolyzer_type == "ALK":
            H2_cell[ii] = electrolyzer.cell_H2_production_rate(I_stack, electrolyzer.T_stack)
        else:
            H2_cell[ii] = electrolyzer.cell_H2_production_rate(I_stack)

    H2_stack = H2_cell*electrolyzer.n_cells
    Stack_Power_kW = current_range*V_cell*electrolyzer.n_cells/1e3
    power_usage_kWh = current_range*V_cell/1e3
    efficiency_kWh_pr_kg = power_usage_kWh/H2_cell
    keys = ["I [A]","V_cell [V]","Stack Power [kW]","Stack H2 Production [kg/hr]","Efficiency [kWh/kg]","Load %"]
    vals = [current_range,V_cell,Stack_Power_kW,H2_stack,efficiency_kWh_pr_kg,100*Stack_Power_kW/electrolyzer.stack_rating_kW]
    df = pd.DataFrame(dict(zip(keys,vals)))
    df.to_csv(filepath)
    print("Saved efficiency curve to: {}".format(filepath))
    return df