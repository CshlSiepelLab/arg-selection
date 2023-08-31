# **S**election **I**nference using the **A**ncestral recombination graph (SIA)
<img src="sia_logo.png" alt="drawing" width="100"/>

This GitHub repo hosts the code for running the SIA method to infer positive selection from genomic polymorphism data as described [here](https://doi.org/10.1093/molbev/msab332).

__N.b. Since the introduction of the domain-adaptive version of SIA $-$ dadaSIA, we recommend users to adopt [this improved method](https://github.com/ziyimo/popgen-dom-adapt) which addresses the issue of simulation mis-specification.__

The SIA framework consists of the following 3 major steps:

## 1 Training data simulation

SIA is a supervised machine learning method that relies on simulated training data. The simulations should be tailored to the population of interest. You can use any simulator of your choice that is able to simulate selection (such as [SLiM](https://messerlab.org/slim/) and [discoal](https://github.com/kr-colab/discoal)). In the `sims` directory, we provide examples of SLiM simulation scripts. See the [documentation](https://github.com/kr-colab/discoal/blob/master/discoaldoc.pdf) for how to run discoal on the command line.

## 2 ARG inference and feature extraction

The SIA deep learning model uses features of the Ancestral Recombination Graph (ARG) for selection inference. Therefore, to train/apply the SIA model, we need to first reconstruct ARGs from polymorphism data.
There are a few methods to choose from for ARG inference, including [ARGweaver](https://github.com/mjhubisz/argweaver), [Relate](https://myersgroup.github.io/relate/) and [tsinfer/tsdate](https://github.com/tskit-dev/tsinfer). We chose to use Relate in our analyses to balance the [tradeoff](https://doi.org/10.1101/2021.11.15.468686) between accuracy and scalability, but the user is free to pick their own ARG inference method.

___There are two versions of ARG features for the SIA model:___

### 2.1 Lineage counts at discretized time points

This is the feature encoding proposed in the [original SIA publication](https://doi.org/10.1093/molbev/msab332) for the base version of SIA. `util.py` contains function `ARG2feature()` for extracting features from the inferred ARG and requires the [DendroPy](https://dendropy.org/) package as a dependency. Import the function by
```python
from util import ARG2feature
```
 and call the function as follows:

```python
ARG2feature(intvls, trees, var_ppos, VOI_gt, no_ft, time_pts, taxa)
```

| Input | Data type | Description |
| ----- | --------- | ----------- |
| `intvls` | numpy array of shape `(k, 2)` and type `int` where `k` is the # of genealogies | Each row encodes the starting (inclusive) and ending (exclusive) position of a genealogy |
| `trees` | list of length `k` | List of newick strings of genealogies |
| `var_ppos` | `int` | Physical position of focal variant (i.e. putative sweep site) |
| `VOI_gt` | numpy array of shape `(m,)` where `m` is the # of haplotypes | Binary array of the genotype at the focal variant, `0` for ancestral allele and `1` for derived allele |
| `no_ft` | `int` | # of flanking trees for feature extraction, we used `2` in analyses presented in the manuscript |
| `time_pts` | numpy array of shape `(n,)` | Array of discretized time points, this can be generated manually or by calling the `discretizeT()` function in `util.py` |
| `taxa` | `int` | The starting index of taxa IDs |

**Output:**
| Data type | Description |
| --------- | ----------- |
| numpy array of shape `(2*no_ft+1, n)` and type `int` | This is the feature matrix of the ARG centered on the variant of interest, see **ARG Feature Extraction** section in [manuscript](https://doi.org/10.1093/molbev/msab332) for details |

### 2.2 Full encoding of genealogies using stacked matrices

We have since introduced a **d**omain-**ada**ptive version of SIA (dadaSIA). Details of the improvements can be found in the [manuscript](https://doi.org/10.1101/2023.03.01.529396). dadaSIA uses a new feature input that is fully representative of all information in the genealogies. Although we recommend users to adopt the new dadaSIA method for selection inference, the improved feature encoding can be used in conjunction with the base version of SIA. This stacked-matrices encoding scheme is implemented via the `encode()` function in `fea_encoding_v2.py`, whose usage is described [here](https://github.com/ziyimo/popgen-dom-adapt#1-domain-adaptive-sia-dadasia).

## 3 Deep learning

The original SIA model uses a stacked LSTM architecture and we provide scripts for training the model in the `DL_models` directory.

* `SIA_swp_class.py` implements the classification model to identify sweeps.
* `SIA_sc_reg.py` implements the regression model to infer selection coefficients.

In order to utilize the new feature encoding with stacked matrices described in [2.2](#22-full-encoding-of-genealogies-using-stacked-matrices), the neural network architecture needs to be modified to use convolutional layers (as described [here](https://doi.org/10.1101/2023.03.01.529396)). These CNN models are implemented in `SIA_swp_class_CNN.py` and `SIA_sc_reg_CNN.py`.

<br>
References:

> Hejase HA, Mo Z, Campagna L, Siepel A. 2021. A Deep-Learning Approach for Inference of Selective Sweeps from the Ancestral Recombination Graph. _Molecular Biology and Evolution_:msab332.

> Mo Z, and Siepel A. 2023. Domain-adaptive neural networks improve supervised machine learning based on simulated population genetic data. _bioRxiv_:10.1101/2023.03.01.529396.