# Bootstrap RMSF
General recipe to compute bootstraped RMSF from MD simulations.

## Introduction

## Documentation

### Prerequisites

You need, at a minimum (requirements.txt):

* Python 2.7
* NumPy
* Pandas
* Matplotlib
* MDAnalysis

### Bootstrap RMSF

```python
# Modules to import at the beginning
from MDAnalysis import Universe

from bootstrap_rmsf import Bootstrap_RMSF
from bootstrap_rmsf import plot_rmsf
```

1. The first step will be to open the trajectories and select the group of atoms for which we want to compute the bootstraped RMSF. Of course, we assume that you are providing trajectories where the protein was aligned on a reference structure. Here in this example, only the heavy atoms that belong to the backbone are selected from the trajectory.

```python
# You can open one trajectory or a list of trajectories
u = Universe('protein.psf', ['trajectory.dcd'])
atomgroup = u.select_atoms(selection='backbone')
```

2. The next step consists to compute the bootstrap mean RMSF for each residue selected. In order to calculate the bootstrap RMSF, *n_sample* conformations will be selected randomly with replacement from the trajectories. The random selection is then repeated *n_iteration* times, giving us a mean RMSF and standard deviation of the mean for each residue selected. Here in this example, 10,000 conformations will be selected randomly and repeated 1,000 times. As result, you will obtain a Pandas DataFrame with the following columns: resid, resname, segid, the mean RMSF and the standard deviation of the mean.

```python
# Do a bootstrap
b = Bootstrap_RMSF(n_sample=10000, n_iteration=1000)
rmsf = b.run(atomgroup)
print rmsf
```

```bash
     resid resname segid      mean       std
1        1     PHE     A  2.984166  0.029518
3        2     TYR     A  2.336458  0.024672
..     ...     ...   ...       ...       ...
458    264     CYS     B  3.072247  0.046174
459    265     PRO     B  3.495771  0.049302
```

3. At the end, you can either save the results in a CSV file and/or directly plot the mean RMSF with the standard deviation associated along the sequence. In that particular exemple, only the RMSF of residues that belong to the segid A will be plotted.

```python
# Save results in a CSV for later use
b.save('rmsf_backbone.csv')
# ... or/and plot the RMSF
plot_rmsf('bootstrap_rmsf_A_CA.png', rmsf[rmsf['segid'] == 'A'], ymax=4)
```

## License
MIT
