import yaml
from importlib import resources


class Config:
    def __init__(self):
        info_file = resources.path('pattools', 'INFO.yaml')
        with open(str(info_file), 'r') as file:
            data = yaml.safe_load(file)
            self.VERSION = data.get('VERSION', '')


CONFIG = Config()
