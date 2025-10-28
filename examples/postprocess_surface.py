#!/bin/env python3
################################################################################
import argparse
import os
from datetime import datetime
from typing import Any

import numpy as np
from h5py import File

from surface import Surfaces
from surface import surface_func as sf

################################################################################

parser = argparse.ArgumentParser(
    description="Postprocessing for surface outputs"
)
parser.add_argument("simpath", type=str, help="Path to simulation")
parser.add_argument("-b", "--batchtools",  action="store_true",
                    help="assume batchtools subdirectory structure (output-0000...)")
parser.add_argument("-o", "--outputpath", default=None,
                    help="Directory to output to")
parser.add_argument("-s", "--isurf", default=1, type=int,
                    help="Index of surface to use")
parser.add_argument("-r", "--irad", default=0, type=int,
                    help="Index of radius to use")
parser.add_argument("-c", "--criteria", nargs="*", default=["bernoulli_out"], type=str,
                    help=("Ejection criteria to compute."
                          " Choices: geodesic, bernoulli, "
                          "geodesic_out, bernoulli_out, none, none_out"))
parser.add_argument("-g", "--histograms",  nargs="*", default=[], type=str,
                    help=("Ejecta histograms to calculate. "
                          " Choices: vinf or any dataset name in the files"))
parser.add_argument("-w", "--weighted_averages",  nargs="*", default=[], type=str,
                    help=("Weighted average time series to calculate. "
                          " Choices: vinf, th, ph, tau or any dataset name in the files"))
parser.add_argument("-m", "--mass_ejection", action="store_true",
                    help="Calculate masse ejection rate")
parser.add_argument("-l", "--nu_luminosities", action="store_true",
                    help="Calculate neutrino luminosities")
parser.add_argument("-e", "--nu_energies", action="store_true",
                    help="Calculate neutrino energies")
parser.add_argument("-t", "--backtrack_temp", type=float, default=5.0,
                    help="Temperature for backtracking in GK. Default=8")
parser.add_argument("-y", "--ejecta", action="store_true",
                    help="Create ejecta.h5 for use with tabulated yields")
parser.add_argument("-E", "--eos", default=None, type=str,
                    help="EOS table in pycompose format")
parser.add_argument("-n", "--ncpu", default=1, type=int,
                    help="Number of cores to use")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Print progress")
parser.add_argument("--t_min", type=float, default=None,
                    help="Minimum time taken into account in analysis")
parser.add_argument("--t_max", type=float, default=None,
                    help="Maximum time taken into account in analysis")

args = parser.parse_args()

if args.outputpath is None:
    today = datetime.today().strftime("%Y-%m-%d")
    args.outputpath = f"{args.simpath}/surface_analysis_s{args.isurf}_r{args.irad}_{today}"

os.makedirs(args.outputpath, exist_ok=True)


################################################################################

if args.batchtools:
    paths = [f"{args.simpath}/{d}" for d in os.listdir(args.simpath)
             if (os.path.isdir(f"{args.simpath}/{d}")
                 and d.startswith("output-"))]
else:
    paths = [args.simpath]

temp_bt = args.backtrack_temp / 11.60452

s = Surfaces(
    paths,
    args.isurf,
    args.irad,
    n_cpu=args.ncpu,
    verbose=args.verbose,
    eos_path=args.eos,
    t_min=args.t_min,
    t_max=args.t_max,
    )
dt = s.times[1] - s.times[0]

################################################################################

bins = {
    "passive_scalars.r_0": np.linspace(0, 0.65, 66),
    "vinf": np.linspace(0, 1, 101),
    "hydro.aux.s": np.linspace(0, 250, 101),
    "hydro.aux.T": np.linspace(0, 1, 101),
    "tau": np.geomspace(20.3, 20300., 101),
    "tau_b": np.geomspace(20.3, 20300., 101),
    "ph": np.linspace(0, 2*np.pi, len(s.aux["ph"])+1),
    "th": np.linspace(0, np.pi, len(s.aux["th"])+1),
    }

msol_to_ms = 0.004925490948309319

# copied from outflowed.cc
def make_bin_from_centers(c: np.ndarray) -> np.ndarray:
    assert c.size > 1, "Need at least two centers"
    b = np.empty(c.size + 1, dtype=c.dtype)
    delta = c[1] - c[0]
    b[0] = c[0] - 0.5 * delta
    for i in range(1, c.size):
        b[i] = c[i-1] + 0.5 * delta
        delta = c[i] - c[i-1]
    b[-1] = c[-1] + 0.5 * delta
    return b

# copied from runs-thc-ba:analysis/ejecta/skynet/grid.h5
ej_ye = np.array((0.01, 0.04, 0.07, 0.1, 0.13, 0.16, 0.19, 0.22, 0.25,
    0.29, 0.32, 0.35, 0.38, 0.41, 0.44, 0.47, 0.5))
ej_s = np.array((1, 1.3, 1.8, 2.4, 3.2, 4.2, 5.6, 7.5, 10, 13, 18, 24,
    32, 42, 56, 75, 100))
ej_tau = np.array((0.1, 0.17, 0.29, 0.49, 0.84, 1.4, 2.4, 4.2, 7.1, 12,
    21, 35, 59, 100, 170, 290, 500))

ejecta_bins = tuple(make_bin_from_centers(ej)
                    for ej in (ej_ye, ej_s, ej_tau/msol_to_ms))

################################################################################

def _check_extra(f: str, crit: str) -> str | sf.SurfaceFunc:
        if f == "vinf":
            if "bernoulli_min" in crit:
                return sf.vinf_min(s.eos)
            elif "bernoulli_eos" in crit:
                return sf.vinf_eos(s.eos)
            elif "bernoulli" in crit:
                return sf.vinf["bernoulli"]
            return sf.vinf["geodesic"]
        if f == "tau":
            return sf.tau
        if f == "tau_b":
            return sf.tau_b(temp_bt, s.eos)
        return f

surf_funcs = {}

ut = {
    "none": None,
    "bernoulli": "hydro.aux.hu_t",
    "bernoulli_min": sf.hut(s.eos),
    "bernoulli_eos": sf.hut_eos(s.eos),
    "geodesic": "hydro.aux.u_t",
    }

for crit in args.criteria:
    out = crit.endswith("_out")
    if out:
        cr = ut[crit[:-4]]
    else:
        cr = ut[crit]
    if args.mass_ejection:
        surf_funcs[f"sc_mdot_{crit}"] = sf.mdot(crit=cr, out=out)
    for f in args.weighted_averages:
        _f = _check_extra(f, crit)
        surf_funcs[f"sc_mdot_{f}_{crit}"] = sf.mdot(weight=_f, crit=cr, out=out)
    for f in args.histograms:
        _f = _check_extra(f, crit)
        surf_funcs[f"hist_{f}_{crit}"] = sf.mass_histogram(_f, out=out, crit=cr, bins=bins[f])
    if args.ejecta:
        surf_funcs[f"ejecta_{crit}"] = sf.ejecta(temp_bt, bins=ejecta_bins, out=out, crit=cr, eos=s.eos)

if args.nu_luminosities:
    for inu in range(3):
        surf_funcs[f"sc_nu{inu}_lum"] = sf.nu_lum[inu]

if args.nu_energies:
    for inu in range(3):
        surf_funcs[f"sc_nu{inu}_en"] = sf.nu_e[inu]

all_funcs = sf.get_many(surf_funcs.values())
all_funcs.name = "surface reductions"

data: dict[str, Any] = {k: [] for k in surf_funcs}
for raw_data in s.process_h5_parallel((all_funcs,), ordered=True):
    for k, d in zip(surf_funcs.keys(), raw_data[0]):
        data[k].append(d)

for c in args.criteria:
    if args.mass_ejection:
        data[f"sc_mej_{c}"] = np.cumsum(data[f"sc_mdot_{c}"])*dt

    for h in args.histograms:
        hist, bin_edges = zip(*data[f"hist_{h}_{c}"])
        data[f"hist_{h}_{c}_cum"] = sum(hist), bin_edges[0]

    if args.ejecta:
        hist, bin_edges = zip(*data[f"ejecta_{c}"])
        data[f"ejecta_{c}_cum"] = sum(hist), bin_edges[0]


################################################################################

scalars = tuple(k[3:] for k in data if  k.startswith("sc_"))

for sc in scalars:
    data[f"sc_{sc}"] = np.array(data[f"sc_{sc}"])

################################################################################

scalar_headers = [sc.replace("bernoulli", "b").replace("geodesic", "g")
                  for sc in scalars]
header = ("{:<24s} "*(len(scalars)+1)).format("time", *scalar_headers)
sdata = np.stack([s.times] + [data[f"sc_{sc}"] for sc in scalars], axis=1)
np.savetxt(f"{args.outputpath}/scalars.txt", sdata, fmt="%24.16e", header=header)

################################################################################

hists = tuple(k[5:] for k in data if k.startswith("hist_") and not k.endswith("_cum"))
for h in hists:
    hist, bin_edges = data[f"hist_{h}_cum"]
    with open(f"{args.outputpath}/hist_{h}.txt", "w") as hf:
        hf.writelines((" ".join(bin_edges.astype(str)) + "\n",
                       " ".join(hist.astype(str)) + "\n"))

################################################################################

if args.ejecta:
  for crit in args.criteria:
    try:
        mass, _ = data[f"ejecta_{crit}_cum"]
    except ValueError:
        print(len(data[f"ejecta_{crit}"]))
        raise
    with File(f"{args.outputpath}/ejecta_{crit}.h5", "w") as hf:
        hf['mass'] = mass
        hf['Ye'] = ej_ye
        hf['entropy'] = ej_s
        hf['tau'] = ej_tau
