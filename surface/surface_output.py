################################################################################

import os
import numpy as np

################################################################################

def _straighten(data: dict) -> dict:
    if 'iter' in data:
        dsort = data['iter']
    elif 'time' in data:
        dsort = data['time']
    else:
        return data
    _, isort = np.unique(dsort, return_index=True)
    return {k: np.atleast_1d(dd)[isort] for k, dd in data.items()}

def _read_ascii(path: str) -> dict:
    with open(path, 'r') as f:
        line = "#"
        header = ""
        while line.startswith("#"):
            header = line
            line = f.readline()
    keys = [hh for hh in header.split()[1:]]

    data = np.loadtxt(path, skiprows=1, unpack=True)
    return _straighten(dict(zip(keys, data)))

def _read_hist(path: str) -> tuple[np.ndarray, np.ndarray]:
    with open(path, 'r') as hf:
        bins = np.array(hf.readline().split()).astype(float)
        hst = np.array(hf.readline().split()).astype(float)
    return hst, bins

################################################################################

class SurfaceOutput:

    def __init__(self, path: str):
        self.path = path
        files = os.listdir(path)
        if "scalars.txt" in files:
            self.scalars = _read_ascii(f"{path}/scalars.txt")
        else:
            self.scalars = {}
        self.hist = {hist[5:-4]: _read_hist(f"{path}/{hist}")
                     for hist in files if hist.startswith("hist_")}
