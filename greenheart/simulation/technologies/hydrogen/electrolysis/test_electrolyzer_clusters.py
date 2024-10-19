from greenheart.simulation.technologies.hydrogen.electrolysis.ALK_electrolyzer_clusters import ALK_Clusters
from greenheart.simulation.technologies.hydrogen.electrolysis.PEM_sync_electrolyzer_clusters import PEM_Clusters

from greenheart.simulation.technologies.hydrogen.electrolysis.ALK_electrolyzer_tools import get_efficiency_curve as alk_curve
from greenheart.simulation.technologies.hydrogen.electrolysis.PEM_electrolyzer_tools import get_efficiency_curve as pem_curve


if __name__ == "__main__":
    alk = ALK_Clusters(cluster_size_mw=1, plant_life=30)
    df = alk_curve(alk, file_desc = "July2024")
    
    print(alk)
    print(df)

    # pem = PEM_Clusters(cluster_size_mw=1, plant_life=30)
    # # df = pem_curve(pem, file_desc = "July2024")
    # for item in dir(pem):
    #     print(item)
    # print(pem)
    # # print(df)