from greenheart.simulation.technologies.hydrogen.electrolysis.alkaline_control_operation import AlkalineSupervisor
import numpy as np


def run_alkaline_physics():
    """1 liner desc

        longer desc:cite:`jvm-jensen1983note`

        Args:
            input_signal_type (str): options are "power" or "h2"
            degradation_application (str): options are "power" or "h2"
            control_strategy (str): options are "even_split_load"

            alkaline_config (dict): _description_
            V_init (_type_): _description_
            V_deg (_type_): _description_

        Returns:
            _type_: _description_
        
        """

    num_clusters = int(np.ceil(electrolyzer_size_MW/cluster_size_MW))
    sup = AlkalineSupervisor(electrolyzer_size_MW,cluster_size_MW,plant_life)
    clusters = sup.create_clusters(num_clusters,cluster_size_MW,alk_config)
    input_signal = "power"
    hydrogen_signal_for_clusters

