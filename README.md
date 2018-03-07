# Bootstrap RMSF
General recipe to compute bootstraped RMSF from MD simulations.

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
```

1. The first step will be to open the trajectories and select the group of atoms for which we want to compute the bootstraped RMSF. Of course, we assume that you are providing trajectories where the protein was aligned on a reference structure. Here in this example, only the heavy atoms that belong to the backbone are selected from the trajectory.

```python
# You can open one trajectory or a list of trajectories
u = Universe('protein.psf', ['trajectory.dcd'])
atomgroup = u.select_atoms(selection='backbone')
```

2. The next step consists to compute the bootstrap RMSF for each residue selected. In order to calculate the bootstrap RMSF, *n_sample* conformations will be selected randomly with replacement from the trajectories. Normally, the number of selected conformations is equal to the total number of conformations. The random selection is then repeated *n_iteration* times, giving us the (pivot) confidence interval around the RMSF for each residue selected. Here in this example, all conformations (```n_sample=-1```) will be selected randomly with replacement and this process is repeated 1,000 times. We estimate the 95% confidence interval. As result, you will obtain a Pandas DataFrame with the following columns: resid, resname, segid, the RMSF and the confidence interval. The results can be saved in a CSV file using the *save* function.

```python
# Do a bootstrap
b = Bootstrap_RMSF(n_iteration=1000, alpha=0.05, n_sample=-1)
rmsf = b.run(atomgroup)

print rmsf

# Save results in a CSV for later use
b.save('rmsf_backbone.csv')
```

```bash
     resid resname segid      rmsf       low      high
1        1     PHE     A  1.634058  1.534535  1.730403
3        2     TYR     A  1.139555  1.070607  1.223800
4        3     CYS     A  0.857350  0.798316  0.920177
..     ...     ...   ...       ...       ...       ...
458    263     GLY     B  1.403246  1.189719  1.609360
459    264     CYS     B  1.391054  1.209021  1.579265
460    265     PRO     B  2.013367  1.837854  2.171748

```

3. (Optional) At the end, you can also plot directly the RMSF along the sequence. In that particular exemple, only the RMSF of residues that belong to the segid A will be plotted. Also, resids will be renumbered with the first residu starting at 507 instead of 1.

```python
from bootstrap_rmsf import plot_rmsf

plot_rmsf('bootstrap_rmsf_A_CA.png', rmsf[rmsf['segid'] == 'A'], ymax=4, start_resid=507)
```

## References
1. [Bootstrap methods: Another look at the jackknife (original paper)](https://projecteuclid.org/download/pdf_1/euclid.aos/1176344552)
2. [The Bootstrap](http://www.stat.cmu.edu/~cshalizi/402/lectures/08-bootstrap/lecture-08.pdf)
3. [Bootstrap confidence intervals](https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf) 

## To-do list
- [ ] Python 3
- [ ] Multiprocessing
- [ ] Welch t-test

## License
MIT
