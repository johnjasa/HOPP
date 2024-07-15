from greenheart.simulation.technologies.hydrogen.electrolysis.alkaline_control_operation import AlkalineSupervisor
import numpy as np

def combine_results_across_clusters(res):
    pass
def run_alkaline_physics(input_signal,input_signal_type,electrolyzer_size_MW,electrolyzer_config):
    """1 liner desc

        longer desc:cite:`jvm-jensen1983note`

        Args:
            input_signal (np.NDarray or float): options are. signal must correspond to input_signal_type
                - hydrogen demand profile in kg/hr (float or array) or 
                - input power in kW
            input_signal_type (str): options are "power" or "h2"
            degradation_application (str): options are "power" or "h2"
            control_strategy (str): options are "even_split"

            alkaline_config (dict): _description_
            V_init (_type_): _description_
            V_deg (_type_): _description_

        Returns:
            _type_: _description_
        
        """
    alk_config = electrolyzer_config["cluster_config"]
    control_config = electrolyzer_config["operation_config"]
    control_strategy = "even_split"
    cluster_size_MW = alk_config["cluster_size_mw"]
    plant_life = alk_config["plant_life"]
    num_clusters = int(np.ceil(electrolyzer_size_MW/cluster_size_MW))
    sup = AlkalineSupervisor(electrolyzer_size_MW,cluster_size_MW,input_signal_type,control_strategy,plant_life)
    
    clusters = sup.create_clusters(num_clusters,cluster_size_MW,alk_config)
    res, power_consumption_total,hydrogen_production_total = sup.run(clusters,input_signal)

    

