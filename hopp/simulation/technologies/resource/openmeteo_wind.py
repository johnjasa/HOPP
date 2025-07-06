import os
from pathlib import Path
from typing import Union, Optional, List
import urllib.parse

from attrs import define, field

from hopp.utilities.keys import get_developer_nrel_gov_key, get_developer_nrel_gov_email
from hopp.utilities.validators import range_val
from hopp.simulation.technologies.resource.resource import Resource
from hopp import ROOT_DIR
from hopp.tools.resource.pysam_wind_tools import combine_wind_files

import matplotlib.pyplot as plt
import numpy as np
import openmeteo_requests
import pandas as pd
import requests_cache
from retry_requests import retry


@define
class MeteoWindData(Resource):    
    """
    """
    
    lat: float = field()
    lon: float = field()
    #: year for resource data. Must be between 1950 and 2025
    year: int = field(validator=range_val(1950, 2025))

    #: the hub-height for wind resource data (meters)
    hub_height_meters: float = field(validator=range_val(10.0, 200.0))
    
    # OPTIONAL INPUTS
    path_resource: Optional[Union[str, Path]] = field(default = ROOT_DIR / "simulation" / "resource_files")
    filename: Optional[Union[str, Path]] = field(default = None)
    use_api: Optional[bool] = field(default = False)
    resource_data: Optional[dict] = field(default = None)

    #: dictionary of heights and filenames to download from Wind Toolkit
    file_resource_heights: dict = field(default = None)

    # NOT INPUTS
    allowed_hub_height_meters: List[int] = [10, 80, 120, 180]
    

    def __attrs_post_init__(self):
        super().__init__(self.lat, self.lon, self.year)

        if self.interval != "60":
            raise ValueError("Only hourly (interval=60 mins) wind data is supported for METEO currently.")

        # if resource_data is input as a dictionary then set_data   
        if isinstance(self.resource_data,dict):
            self.data = self.resource_data
            return 
        
        # if resource_data is not provided, download or load resource data
        if isinstance(self.path_resource,str):
            self.path_resource = Path(self.path_resource).resolve()
        if self.path_resource.parts[-1]!="wind":
            self.path_resource = self.path_resource / 'wind'

        if self.filename is None:
            self.calculate_heights_to_download()

        self.check_download_dir()

        if not os.path.isfile(self.filename) or self.use_api:
            responses = self.download_resource()
        
        self.format_data(responses)
        
    def calculate_heights_to_download(self):
        """
        Given the system hub height, and the available hubheights from METEO Data,
        determine which heights to download to bracket the hub height
        """
        hub_height_meters = self.hub_height_meters

        # evaluate hub height, determine what heights to download
        heights = [hub_height_meters]
        if hub_height_meters not in self.allowed_hub_height_meters:
            height_low = self.allowed_hub_height_meters[0]
            height_high = self.allowed_hub_height_meters[-1]
            for h in self.allowed_hub_height_meters:
                if h < hub_height_meters:
                    height_low = h
                elif h > hub_height_meters:
                    height_high = h
                    break
            heights[0] = height_low
            heights.append(height_high)

        filename_base = f"{self.latitude}_{self.longitude}_METEO_{self.year}_{self.interval}min"
        file_resource_full = filename_base
        file_resource_heights = dict()

        for h in heights:
            h_int = int(h)
            file_resource_heights[h_int] = self.path_resource/(filename_base + f'_{h_int}m.csv')
            file_resource_full += f'_{h_int}m'
        file_resource_full += ".csv"

        self.file_resource_heights = file_resource_heights
        self.filename = self.path_resource / file_resource_full

    def update_height(self, hub_height_meters):
        self.hub_height_meters = hub_height_meters
        self.calculate_heights_to_download()

    def download_resource(self):
        """
        Downloads the wind data from the METEO dataset using an API call
        """
        success = False

        base_attributes = ["temperature_2m", "wind_speed_10m", "wind_direction"]
        attributes = ["pressure_msl"]
        for height, f in self.file_resource_heights.items():
            attributes += [f"{a}_{height}m" for a in base_attributes]
        
        # Setup the Open-Meteo API client with cache and retry on error
        cache_session = requests_cache.CachedSession(".cache", expire_after=3600)
        retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
        openmeteo = openmeteo_requests.Client(session=retry_session)

        # Make sure all required weather variables are listed here
        # The order of variables in hourly or daily is important to assign them correctly below
        url = "https://historical-forecast-api.open-meteo.com/v1/forecast"
        params = {
            "latitude": self.latitude,
            "longitude": self.longitude,
            "start_date": f"{self.year}-01-01",
            "end_date": f"{self.year}-12-31",
            "hourly": attributes,
            "wind_speed_unit": "ms",
            "temperature_unit": "celsius",
        }

        # PAUL NOTES: GET SSL ERRORS UNLESS I INCLUDE verify=False
        # BUT UNDERSTAND THIS IS NOT A GOOD PRACTICE
        responses = openmeteo.weather_api(url, params=params, verify=False)

        if not responses:
            raise ValueError('Unable to download wind data')

        return responses

    def format_data(self):
        """
        Format as 'wind_resource_data' dictionary for use in PySAM.
        """
        self.data = self.filename

    @Resource.data.setter
    def data(self, data_info):
        """
        Sets the wind resource data to a dictionary in SAM Wind format (see Pysam.ResourceTools.SRW_to_wind_data)
        """
        if isinstance(data_info,dict):
            self._data = data_info
        if isinstance(data_info,(str, Path)):
            resource_heights = [k for k in self.file_resource_heights.keys()]
            self._data = combine_wind_files(str(data_info),resource_heights)