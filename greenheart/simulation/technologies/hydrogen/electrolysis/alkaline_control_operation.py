from greenheart.simulation.technologies.hydrogen.electrolysis.ALK_electrolyzer_clusters import ALK_Clusters as alk
from greenheart.simulation.technologies.hydrogen.electrolysis.alkaline_LTA import alkaline_LTA
import numpy as np
import pandas as pd
class AlkalineSupervisor:
    def __init__(self,electrolyzer_size_MW,cluster_size_MW,plant_life):
        self.plant_life = plant_life
        # if alkaline_config is None:
        #     alk_kwargs = {}
        # else:
        #     alk_kwargs = alkaline_config
        num_clusters = int(np.ceil(electrolyzer_size_MW/cluster_size_MW))
        self.electrolyzer_size_MW = num_clusters*cluster_size_MW

        # clusters = self.create_clusters(num_clusters,cluster_size_MW,alk_kwargs)

    def run(self,clusters):
        for ci, cluster in enumerate(clusters):
            clusters[ci].run_cluster_variable_power(power_to_clusters[ci])
            clusters[ci].run_cluster_hydrogen_demand(H2_required_cluster_kg[ci])
            clusters[ci].simulation_results
            clusters[ci].timeseries_results
            clusters[ci].BOL_design_info
            if clusters[ci].run_LTA:
                clusters[ci].LTA_results_annual


        # power_to_clusters = self.even_split_power(input_power_kW,clusters)
        # h2_df_ts = pd.DataFrame()
        # h2_df_tot = pd.DataFrame()
        # col_names = []
        # for ci, cluster in enumerate(clusters):
        #     cl_name = "Cluster #{}".format(ci)
        #     col_names.append(cl_name)
        #     h2_ts, h2_tot = clusters[ci].run(power_to_clusters[ci])
        #     # h2_dict_ts['Cluster #{}'.format(ci)] = h2_ts

        #     h2_ts_temp = pd.Series(h2_ts, name=cl_name)
        #     h2_tot_temp = pd.Series(h2_tot, name=cl_name)
        #     if len(h2_df_tot) == 0:
        #         h2_df_tot = pd.concat(
        #             [h2_df_tot, h2_tot_temp], axis=0, ignore_index=False
        #         )
        #         h2_df_tot.columns = col_names

        #         h2_df_ts = pd.concat([h2_df_ts, h2_ts_temp], axis=0, ignore_index=False)
        #         h2_df_ts.columns = col_names
        #     else:
        #         # h2_df_ts = h2_df_ts.join(h2_ts_temp)
        #         h2_df_tot = h2_df_tot.join(h2_tot_temp)
        #         h2_df_tot.columns = col_names

        #         h2_df_ts = h2_df_ts.join(h2_ts_temp)
        #         h2_df_ts.columns = col_names
        # pass

    def even_split_power_sequential(self):
        pass

    def even_split_hydrogen(self,hydrogen_demand_profile_kg,clusters):
        #NOTE: assumes all clusters are the same rating and capcaity
        i=0
        num_clusters_on = np.floor(hydrogen_demand_profile_kg / clusters[i].cluster_min_h2_kg)
        num_clusters_on = np.where(
            num_clusters_on > len(clusters), len(clusters), num_clusters_on
        )
        hydrogen_per_cluster = [
            hydrogen_demand_profile_kg[ti] / num_clusters_on[ti]
            if num_clusters_on[ti] > 0
            else 0
            for ti, pwr in enumerate(hydrogen_demand_profile_kg)
        ]
        hydrogen_to_active_clusters = np.array(hydrogen_per_cluster)
        hydrogen_signal_for_clusters = np.zeros((hydrogen_demand_profile_kg, len(clusters)))
        for i, cluster_power in enumerate(hydrogen_to_active_clusters):
            clusters_off = len(clusters) - int(num_clusters_on[i])
            no_hydrogen = np.zeros(clusters_off)
            with_hydrogen = cluster_power * np.ones(int(num_clusters_on[i]))
            tot_hydrogen = np.concatenate((with_hydrogen, no_hydrogen))
            hydrogen_signal_for_clusters[i] = tot_hydrogen


        return np.transpose(hydrogen_signal_for_clusters)
    def even_split_power(self,input_power_kW,clusters):
        #NOTE: assumes all clusters are the same rating and capcaity
        i=0
        clusters[i].cluster_rating_kW
        clusters[i].min_cluster_power_kW
        num_clusters_on = np.floor(input_power_kW / clusters[i].min_cluster_power_kW)
        num_clusters_on = np.where(
            num_clusters_on > len(clusters), len(clusters), num_clusters_on
        )
        power_per_cluster = [
            input_power_kW[ti] / num_clusters_on[ti]
            if num_clusters_on[ti] > 0
            else 0
            for ti, pwr in enumerate(input_power_kW)
        ]
        power_to_active_clusters = np.array(power_per_cluster)
        power_signal_for_clusters = np.zeros((input_power_kW, len(clusters)))
        for i, cluster_power in enumerate(power_to_active_clusters):
            clusters_off = len(clusters) - int(num_clusters_on[i])
            no_power = np.zeros(clusters_off)
            with_power = cluster_power * np.ones(int(num_clusters_on[i]))
            tot_power = np.concatenate((with_power, no_power))
            power_signal_for_clusters[i] = tot_power


        return np.transpose(power_signal_for_clusters)

    def even_split_hydrogen_demand_sequential(self):
        pass

    def create_clusters(self,num_clusters:int,cluster_size_mw,alk_kwargs):
        
        clusters = []
        if isinstance(cluster_size_mw,(int,float)):
            for i in num_clusters:
                clusters.append(alk(cluster_size_mw,self.plant_life,**alk_kwargs))
        else:
            for i,cluster_size in enumerate(cluster_size_mw):
                clusters.append(alk(cluster_size))
        return clusters

if __name__ == "__main__":
    num_clusters = 2