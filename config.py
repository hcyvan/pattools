import os
import yaml
from pathlib import Path


class Config:
    def __init__(self, config_path):
        self.config_abs_path = os.path.abspath(config_path)
        # cf = configparser.ConfigParser()
        # cf.read(self.config_abs_path, encoding='utf-8')
        # self.DataDir = cf['DEFAULT'].get("DataDir")
        with open(self.config_abs_path, 'r') as f:
            yaml_data = yaml.safe_load(f)
        self.DataDir = Path(yaml_data['DataDir'])
        self.DataRaw = self.DataDir / 'raw'
    def __str__(self):
        rets = []
        for k, v in vars(self).items():
            if 'Password' in k:
                v = "**********"
            rets.append("{}:{}".format(k, v))
        return "; ".join(rets)


CONFIG = Config(os.path.join(os.path.dirname(__file__), "config.yaml"))
