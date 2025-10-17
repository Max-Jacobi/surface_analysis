# surface_analysis
Postprocessing scripts for surface dumps in GRAthena++

## Installation
Clone and install with `cd surface_analysis; pip install .`
The code requires the [tabulatedEOS](https://github.com/max-jacobi/tabulatedEOS) library which pip will install from github.
If this is not possible (i.e. on clusters), you can clone it and install it by hand (`cd tabulatedEOS; pip install .`).

## Example scripts
Example usage script in [examples](examples/postprocess_surface.py).

To use [yields.py](examples/yields.py), you need the `tabulated_nucsyn.h5` and `solar_r.dat` which e.g. in the runs-thc-ba repo.

``` shell
python ~/repos/surface_analysis/examples/postprocess_surface.py \
    <path/to/sim> \
    -b \ # if batchtools structure
    -v -s1 -r0 -n24 \
    -l -e -m \
    -c bernoulli bernoulli_out geodesic geodesic_out none none_out \ # ejection criteria
    -g vinf hydro.aux.T hydro.aux.s passive_scalars.r_0 tau ph th \ # histograms
    -w vinf hydro.aux.T hydro.aux.s passive_scalars.r_0 tau ph th \ # running weighted averages
    -y -t5 -E <path/to/EOS/file> \ # for nucleosynthesis yield interpolation (might be slowish)
    -o surface_analysis
```
