# surface_analysis
Postprocessing scripts for surface dumps in GRAthena++

## Installation
Clone and install with `cd surface_analysis; pip install .`
The code requires the [tabulatedEOS](https://github.com/max-jacobi/tabulatedEOS) library which pip will install from github.
If this is not possible (i.e. on clusters), you can clone it and install it by hand (`cd tabulatedEOS; pip install .`).

## Example scripts
Example usage script in [examples](examples/postprocess_surface.py).

To use [yields.py](examples/yields.py), you need the `tabulated_nucsyn.h5` and `solar_r.dat` which e.g. in the runs-thc-ba repo.
Full example with all options activated
``` shell
python examples/postprocess_surface.py \
    <path/to/sim> \
    --batchtools \
    --verbose --isurf=1 --irad=0 --ncpu=<number of cpu> \
    --nu_luminosities --nu_energies \
    --mass_ejection \
    --criteria bernoulli bernoulli_out geodesic geodesic_out none none_out \
    --histograms vinf hydro.aux.T hydro.aux.s passive_scalars.r_0 tau ph th \
    --weighted_averages vinf hydro.aux.T hydro.aux.s passive_scalars.r_0 tau ph th \
    --ejecta --eos=<path/to/EOS/file> \
    --t_min=<minimum_time> --t_max=<maximum_time> \
    --outputpath <path/to/output>
```
