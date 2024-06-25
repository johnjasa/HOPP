import math
import numpy as np
import sys
import pandas as pd
from matplotlib import pyplot as plt
import scipy
import rainflow
from scipy import interpolate
from scipy.constants import R, physical_constants, convert_temperature

class ALK_Clusters:
    # num hydrogen molecules transferred per reaction
    z: int = 2 #TODO: change to z_c
    F: float = 96485.34  # Faraday's Constant (C/mol) or [As/mol]
    R: float = 8.314  # Ideal Gas Constant (J/mol/K)

    M_H: float = 1.00784  # molecular weight of Hydrogen [g/mol]
    M_H2: float = 2.016 #[g/mol]
    M_O: float = 15.999  # molecular weight of Oxygen [g/mol]
    M_O2: float = 31.999 #[g/mol]
    M_K: float = 39.0983  # molecular weight of Potassium [g/mol]
    
    lhv: float = 33.33  # lower heating value of H2 [kWh/kg]
    hhv: float = 39.41  # higher heating value of H2 [kWh/kg]
    gibbs: float = 237.24e3  # Gibbs Energy of global reaction (J/mol)
    def __init__(self,cluster_size_mw,plant_life):
        self.dt = 3600 #sec/timestep
        include_degradation_penalty = []
        eol_eff_percent_loss = 10
        uptime_hours_until_eol = 77600

        # OPERATIONAL CONSTRIANTS
        cell_max_current_density = 0.3 #[A/cm^2]
        ramp_rate = [] #percent of rated power per unit time
        
        turndown_ratio = 0.25 
        warm_up_delay = [] #sec to warm up

        # OPERATING CONDITIONS
        self.pressure_operating = 1 # [bar] operating pressure
        T_stack = 60 #check value

        # CLUSTER DESIGN PARAMETERS
        self.n_stacks = []
        # STACK DESIGN PARAMETERS
        self.stack_rating_kW = 1000  # 1 MW
        self.n_cells = []

        # CELL DESIGN PARAMETERS
        self.cell_area = 300 # [cm^2] membrane and electrode area
        self.d_em = 0.125 # [cm] electrode-membrane distance
        self.e_m = 0.05 #membrane thickness - check if used

        self.w_koh = 30 # [wt. %] can range from [25-33]
        # self.electrolyte_concentration_percent = self.w_koh / 100
        

        # CELL DEGRADATION RATES
        self.onoff_deg_rate =  3.0726072607260716e-04 #[V/off-cycle]
        self.rate_fatigue = 1.2820512820512823e-05 #multiply by rf_track
        self.steady_deg_rate = 5.092592592592592e-09 #V/sec

        # INITIALIZATION
        self.m,self.M = self.create_electrolyte()
        self.curve_coeff=self.create_power_current_curve()
# -------------------------------------------- #      
# ----- OPERATIONAL CONSTRAINTS & LOSSES ----- #
# -------------------------------------------- #     
    def cluster_warm_up_losses(self):
        pass
    def check_ramp_rate(self):
        pass
    def run_cluster_variable_power(self):
        pass
    def run_cluster_hydrogen_demand(self):
        pass
    def run_cluster(self,I_stack):
        pass
# ------------------------------------- #      
# ----- CLUSTER - LEVEL EQUATIONS ----- #
# ------------------------------------- #     

    def calc_cluster_status(self,I_stack):
        pass
    def cluster_calc_curtailed_power(self,input_power_kW):
        power_curtailed_kW = np.where(input_power_kW > self.cluster_rating_kW,\
        input_power_kW - self.cluster_rating_kW,0)

        input_power_kW = np.where(input_power_kW >
                        (self.cluster_rating_kW),
                        (self.cluster_rating_kW),
                        input_power_kW)
        pass
# ----------------------------------- #
# ----- STACK - LEVEL EQUATIONS ----- #
# ----------------------------------- #
    def run_stack(self):
        pass
    def create_power_current_curve(self):
        #NOTE: should this be moved to higher level (like AlkalineSupervisor?)
        pass
    def stack_degraded_current(self,I_stack,V_init,V_deg):
        """1 liner desc

        longer desc:cite:`jvm-jensen1983note`

        Args:
            I_stack (_type_): _description_
            V_init (_type_): _description_
            V_deg (_type_): _description_

        Returns:
            _type_: _description_
        
        """
        #current decrease - same power
        # I_in = calc_current((power_input_kW,self.T_C), *self.curve_coeff)
        eff_mult = (V_init + V_deg)/V_init #(1 + eff drop)
        I_deg = I_stack/eff_mult

        return I_deg
    def stack_degraded_power(self,h2_demand_kg,V_init,V_deg):
        #power increase - same current
        pass
# ---------------------------------- #
# ----- CELL - LEVEL EQUATIONS ----- #
# ---------------------------------- #
    def create_electrolyte(self,):
        electrolyte_concentration_percent = self.w_koh / 100
        solution_weight_g = 1000
        density_of_water = 1  # [g/mL] #TODO: could be temperature dependent
        density_of_KOH = 2.12  # [g/mL] #TODO: could be temperature dependent

        self.M_KOH = self.M_O + self.M_H + self.M_K  # [g/mol]
        grams_of_solute = solution_weight_g * (electrolyte_concentration_percent)
        moles_of_solute = grams_of_solute / self.M_KOH  # [mols of solute / solution]
        # solvent is water
        self.M_H2O = 2 * self.M_H + self.M_O  # [g/mol]
        grams_of_solvent = solution_weight_g * (
            1 - electrolyte_concentration_percent
        )
        kg_of_solvent = (1 / 1000) * grams_of_solvent

        volume_of_water = grams_of_solvent / density_of_water  # mL
        volume_of_KOH = grams_of_solute / density_of_KOH  # mL
        volume_of_solution = (volume_of_water + volume_of_KOH) / 1000  # L

        # molality = mol of solute / kg of solvent
        molality = moles_of_solute / kg_of_solvent  # mol/kg
        # molarity is mol/kg = mol of solute / L of solution
        molarity = moles_of_solute / volume_of_solution  # mol/L
        # % solution = amount of solute / amount of solution

        self.m = molality  # NOTE: THIS HAS BEEN VALIDATED
        self.M = molarity  # NOTE: THIS HAS BEEN VALIDATED

        return molality,molarity
    
    
    def cell_bubble_rate_coverage(self, T_stack, I_stack):
        """_summary_

        Args:
            T_stack (_type_): _description_
            I_stack (_type_): _description_

        Returns:
            _type_: _description_
        """
        T_k = convert_temperature([T_stack], "C", "K")[0]
        T_amb = convert_temperature([25], "C", "K")[0]
        J_lim = 30  # [A/cm^2] [Vogt,Balzer 2005]
        j = I_stack / self.cell_area  # [A/cm^2] "nominal current density"
        Pv_H20 = np.exp(
            81.6179 - (7699.68 / T_k) - (10.9 * np.log(T_k)) + (T_k * (9.5891 * 1e-3))
        )
        theta = (
            (self.pressure_operating / (self.pressure_operating - Pv_H20))
            * (-97.25 + 182 * (T_k / T_amb) - 84 * ((T_k / T_amb) ** 2))
            * (j / J_lim) ** (0.3)
        )
        epsilon = (2 / 3) * theta  # bulk bubbling
        return theta, epsilon

    def calc_current_density(self, T_stack, I_stack):
        """_summary_

        Args:
            T_stack (_type_): _description_
            I_stack (_type_): _description_

        Returns:
            _type_: _description_
        """
        theta,epsilon = self.cell_bubble_rate_coverage(T_stack, I_stack)
        A_electrode_eff = self.cell_area * (1 - theta)  # [cm^2]
        j = I_stack / A_electrode_eff  # [A/cm^2]
        return j
    # ---------------------------------- #
    # ----- CELL VOLTAGE EQUATIONS ----- #
    # ---------------------------------- #
    def cell_design(self,T_stack,I_stack):
       
        V_rev = self.cell_reversible_overpotential(T_stack, self.pressure_operating)
        V_act_a, V_act_c = self.cell_activation_overpotential(T_stack, I_stack)
        V_ohm = self.cell_ohmic_overpotential(T_stack, I_stack)
        V_cell = V_rev + V_ohm + V_act_a + V_act_c  # Eqn 4
        return V_cell

    def cell_reversible_overpotential(self, T_stack, P):
        T_K = convert_temperature([T_stack], "C", "K")[0]
        # Eqn 17
        a = (
            (-0.0151 * self.m)
            - (1.6788 * (1e-3) * (self.m**2))
            + (2.2588 * (1e-5) * (self.m**3))
        )
        # Eqn 18 
        b = (
            1
            - ((1.26062 * (1e-3)) * self.m)
            + ((5.6024 * 1e-4) * (self.m**2))
            - ((self.m**3) * (7.8228 * 1e-6))
        )
        # Pv_H20 is vapor pressure of pure water [bar]
        # Eqn 19
        Pv_H20 = np.exp(
            81.6179 - (7699.68 / T_K) - (10.9 * np.log(T_K)) + (T_K * (9.5891 * 1e-3))
        )
        # alpha_h20: water activity of electrolyte solution based on molality [mol/kg]
        # valid for molality ranging from 2-18 mol/kg
        # Eqn 20
        alpha_h20 = np.exp(
            (-0.05192 * self.m)
            + (0.003302 * (self.m**2))
            + (((3.177 * self.m) - (2.131 * self.m**2)) / T_K)
        )
        # Pv_KOH: [bar] vapor pressure of KOH solution
        # Eqn 16
        Pv_KOH = np.exp((2.302 * a) + (b * np.log(Pv_H20)))

        Urev0 = self.cell_Urev0()

        # Eqn 14
        U_rev = Urev0 + ((R * T_K) / (2 * self.F)) * np.log(
            ((P - Pv_KOH) ** 1.5) / alpha_h20
        )
        return U_rev

    def cell_Urev0(self):
        return self.gibbs / (2 * self.F)

    def cell_activation_overpotential(self, T_stack, I_stack):
        # validated against Figure 5 of Reference
        j_eff = self.calc_current_density(T_stack, I_stack)
        ja = j_eff  # [A/cm^2]
        jc = j_eff  # [A/cm^2]
        T_anode = convert_temperature([T_stack], "C", "K")[0]
        T_cathode = convert_temperature([T_stack], "C", "K")[0]
        # Eqn 14 anode charge transfer coeff
        alpha_a = 0.0675 + 0.00095 * T_anode
        # Eqn 15 cathode charge transfer coeff
        alpha_c = 0.1175 + 0.00095 * T_cathode

        # Table 1
        delta_Ga = 41500  # [J/mol*K]
        delta_Gc = 23450  # [J/mol*K]
        jref_0a = 1.34535 * 10 ** (-5)  # [A/cm^2]
        jref_0c = 1.8456 * 10 ** (-3)  # [A/cm^2]
        Tref = convert_temperature([25], "C", "K")[0]
        gamma_a = 1.25  # anode roughness factor
        gamma_c = 1.05  # cathode roughness factor
        # Eqn 16
        j0c = (
            gamma_c
            * jref_0c
            * np.exp((-1 * delta_Gc / R) * ((1 / T_cathode) - (1 / Tref)))
        )
        # Eqn 16
        j0a = (
            gamma_a
            * jref_0a
            * np.exp((-1 * delta_Ga / R) * ((1 / T_anode) - (1 / Tref)))
        )
        # Eqn 13 - Tafel slope for anode
        ba = (R * T_anode) / (self.z * self.F * alpha_a)
        # Eqn 13 - Tafel slope for cathode
        bc = (R * T_anode) / (self.z * self.F * alpha_c)
        # Eqn 11 - anode activation energy
        V_act_a = ba * np.maximum(0, np.log(ja / j0a))
        # Eqn 12 - cathode activation energy
        V_act_c = bc * np.maximum(0, np.log(jc / j0c))

        return V_act_a, V_act_c

    def cell_ohmic_overpotential(self, T_stack, I_stack):
        R_tot = self.cell_total_resistance(T_stack, I_stack)  # Ohms
        V_ohm = I_stack * R_tot  # [V/cell]
        return V_ohm

    def cell_total_resistance(self, T_stack, I_stack):
        R_a, R_c = self.cell_electrode_resistance(T_stack)
        R_electrode = R_a + R_c
        R_ele_bf, R_ele_b = self.cell_electrolyte_resistance(T_stack, I_stack)  # [Ohms]
        R_electrolyte = R_ele_bf + R_ele_b
        R_membrane = self.cell_membrane_resistance(T_stack)  # [Ohms] VERIFIED for Ohm*cm^2
        R_tot = R_electrode + R_electrolyte + R_membrane  # Ohm

        return R_tot
        
    def cell_electrolyte_resistance(self,T_stack, I_stack):
        T_K = convert_temperature([T_stack], "C", "K")[0]
        sigma_bf = (
            -204.1 * self.M
            - 0.28 * self.M**2
            + 0.5332 * self.M * T_K
            + 20720 * (self.M / T_K)
            + 0.1043 * self.M**3
            - 0.00003 * (self.M**2 * T_K**2)
        )
        R_ele_bf = (100 / sigma_bf) * (
            (self.d_em / self.cell_area) + (self.d_em / self.cell_area)
        )
        theta,epsilon = self.cell_bubble_rate_coverage(T_stack, I_stack)
        R_ele_b = R_ele_bf * ((1 / ((1 - epsilon) ** (3 / 2))) - 1)

        return R_ele_bf, R_ele_b  # Ohms

    def cell_membrane_resistance(self, T_stack):
        Rmem = (0.06 + 80 * np.exp(T_stack / 50)) / (
            10000 * self.cell_area
        )  # Equation 36 - Ohms
        return Rmem

    def cell_electrode_resistance(self,T_stack):
        """_summary_

        Args:
            T_stack (_type_): _description_

        Returns:
            Ra (_type_): _description_
            Rc (_type_): _description_
        """
        tref = 25
        temp_coeff = 0.00586  # 1/degC
        # resistivity of 100% dense electrode at tref
        rho_nickle_0 = 6.4 * 10 ** (-6)  # [Ohm*cm]
        # porosity of electrode
        epsilon_Ni = 0.3
        # Eqn 21 - effective resistance of electrode
        rho_nickle_eff = rho_nickle_0 / ((1 - epsilon_Ni) ** 1.5)
        Ra = (
            rho_nickle_eff
            * (self.e_a / self.cell_area)
            * (1 + (temp_coeff * (T_stack - tref)))
        )
        Rc = (
            rho_nickle_eff
            * (self.e_c / self.cell_area)
            * (1 + (temp_coeff * (T_stack - tref)))
        )

        return Ra, Rc  
    
    # -------------------------------------- #
    # ----- CELL DEGRADATION EQUATIONS ----- #
    # -------------------------------------- #
    def cell_degradation(self):
        V_cell = V_cell*cluster_status
        
        V_deg_uptime = self.cell_steady_degradation(V_cell)
        V_deg_onoff = self.cell_onoff_degradation()
        V_fatigue = self.cell_fatigue_degradation()

        V_deg = np.cumsum(V_deg_uptime) + np.cumsum(V_deg_onoff) + V_fatigue
        return V_deg

    def cell_fatigue_degradation(self):
        self.rate_fatigue
        pass

    def cell_steady_degradation(self):
        steady_deg_per_hr=self.dt*self.steady_deg_rate*V_cell*cluster_status
        # cumulative_Vdeg=np.cumsum(steady_deg_per_hr)
        # self.steady_deg_rate
        return steady_deg_per_hr

    def cell_onoff_degradation(self):
        change_stack=np.diff(cluster_status)
        cycle_cnt = np.where(change_stack < 0, -1*change_stack, 0)
        cycle_cnt = np.array([0] + list(cycle_cnt))
        self.off_cycle_cnt = cycle_cnt
        stack_off_deg_per_hr= self.onoff_deg_rate*cycle_cnt
        return stack_off_deg_per_hr
    
    
    # --------------------------------- #
    # ----- CELL OUTPUT EQUATIONS ----- #
    # --------------------------------- #
    def calc_faradaic_efficiency(self,T_stack,I_stack):

        f1 = 250  # [mA^2/cm^4]
        f2 = 0.96  # [-]

        j = self.calc_current_density(T_stack, I_stack)  # [A/cm^2]
        j *= 1000  # [mA/cm^2]
        #Faradaic Efficiency
        eta_F = f2 * (j**2) / (f1 + j**2)
        return eta_F

    def cell_H2_production_rate(self,T_stack,I_stack):
        eta_F = self.calc_faradaic_efficiency(T_stack, I_stack)
        h2_prod_mol = eta_F * I_stack / (2 * self.F) #[mol/sec]
        mfr = self.M_H2 * h2_prod_mol  # [g/sec]
        mfr_H2 = self.dt*mfr / 1e3  # [kg/cell-dt]
        return mfr_H2

    def cell_O2_production_rate(self,T_stack,I_stack):
        eta_F = self.calc_faradaic_efficiency(T_stack, I_stack)
        o2_prod_mol = eta_F * I_stack / (4 * self.F) #[mol/sec]
        mfr = self.M_O2 * o2_prod_mol  # [g/sec]
        mfr_O2 = self.dt*mfr / 1e3  # [kg/cell-dt]
        return mfr_O2
    
    # ------------------------------------ #      
    # ----- ANALYSIS/POST-PROCESSING ----- #
    # ------------------------------------ #     