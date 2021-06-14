# Model-agnostic Feature Importance and Effects with Dependent Features -- A Conditional Subgroup Approach

This repository contains the code, files and data for the JMLR submission "Model-agnostic Feature Importance and Effects with Dependent Features--A Conditional Subgroup Approach"


## Abstract

The interpretation of feature importance in machine learning models is challenging when features are dependent.
Permutation feature importance (PFI) ignores such dependencies, which can cause misleading interpretations due to extrapolation.
A possible remedy is more advanced conditional PFI approaches that enable the assessment of feature importance conditional on all other features. 
Due to this shift in perspective and in order to enable correct interpretations, it is therefore important that the conditioning is transparent and humanly comprehensible.
In this paper, we propose a a new sampling mechanism for the conditional distribution based on permutations in conditional subgroups.
As these subgroups are constructed using decision trees (transformation trees), the conditioning becomes inherently interpretable.
This not only provides a simple and effective estimator of conditional PFI, but also local PFI estimates within the subgroups.
In addition, we apply the conditional subgroups approach to partial dependence plots (PDP), a popular method for describing feature effects that can also suffer from extrapolation when features are dependent and interactions are present in the model.
We show that PFI and PDP based on conditional subgroups often outperform methods such as conditional PFI based on knockoffs, or accumulated local effect plots.
Furthermore, our approach allows for a more fine-grained interpretation of feature effects and importance within the conditional subgroups.


## Reproduce

All experiments are implemented with the R language and the paper is written with LaTeX.

Assuming you have R installed on your system, you can install the package dependencies with R via:

```r
install.packages("devtools")
devtools::install_dev_deps()
```

To reproduce the paper, go to the folder `paper/` and run follwing command in the shell:

```bash
cd paper
make paper
```


## Folder Contents

- `./`: Contains this README and the DESCRIPTION file that specify the R package dependencies
- `./data`: Stores data used in experiments and application.
- `./experiments`: Contains the scripts to produce the figures and results
- `./paper`: Contains the tex files for the paper and the Makefile
- `./R`: Contains custom R functions used in the experiments
- `./results`: Stores intermediate results

## License

&copy; 2020 [Christoph Molnar](https://christophm.github.io/)

The code of this repository is distributed under the MIT license. See
below for details:

```
The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The content of the paper is distributed under a Creative Commons license [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).


