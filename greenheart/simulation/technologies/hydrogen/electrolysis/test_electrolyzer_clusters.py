import unittest
import pandas as pd
from greenheart.simulation.technologies.hydrogen.electrolysis.electrolyzer_base_class import ElectrolyzerCluster

from greenheart.simulation.technologies.hydrogen.electrolysis.ALK_electrolyzer_tools import get_efficiency_curve as alk_curve
from greenheart.simulation.technologies.hydrogen.electrolysis.PEM_electrolyzer_tools import get_efficiency_curve as pem_curve


class TestElectrolyzerCurves(unittest.TestCase):
    
    def test_alk_curve(self):
        # Load the saved ALK truth dataframe
        truth_df = pd.read_pickle('alk_curve.pkl')
        
        # Generate the ALK dataframe
        alk = ElectrolyzerCluster(cluster_size_mw=1, plant_life=30, electrolyzer_type="ALK")
        generated_df = alk_curve(alk, file_desc="July2024")

        pd.testing.assert_frame_equal(generated_df, truth_df)

    def test_pem_curve(self):
        # Load the saved PEM truth dataframe
        truth_df = pd.read_pickle('pem_curve.pkl')
        
        # Generate the PEM dataframe
        pem = ElectrolyzerCluster(cluster_size_mw=1, plant_life=30, electrolyzer_type="PEM")
        generated_df = pem_curve(pem, file_desc="July2024")

        pd.testing.assert_frame_equal(generated_df, truth_df)


if __name__ == "__main__":
    unittest.main()
