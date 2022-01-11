To perform simulations using SLiM:

**1. Make a `.param` file containing demographic parameters for simulations**

The format of the `.param` file is as follows:
```
μ
ρ
gen_0 Ne_0
gen_1 Ne_1
...
```
where `μ` is the **scaled** mutation rate, `ρ` is the **scaled** recombination rate, `gen_i` is the **scaled** generation going **forward** and `Ne_i` is the **scaled** effective population size at `gen_i`.  See `dem_eg.param` for an example.

**2. Neutral and sweep forward simulations**

Use `neutral.slim` and `sweep_treeseq.slim` to run neutral and sweep simulations, respectively. You will need to define the following constants when running the simulations:

| | |
|-|-|
| `paramF` | Name of the `.param` file specifying demographic parameters |
| `outPref` | The file name prefix of the output `.trees` file |

For sweep simulations, these additional constants need to be defined:
| | |
|-|-|
| `selcoef` | The **scaled** selection coefficient |
| `mutgen` | The generation at which the beneficial mutation arises |
| `min_AF` | The minimum derived allele frequency of the beneticial variant at the end of the simulation |
| `max_AF` | The maximum derived allele frequency of the beneticial variant at the end of the simulation |


**3. Recapitation**

Use `recapitation.py` to perform recapitation for the SLiM simulations. Run 
```bash
./recapitation.py
```
for usage.