# Fast non-linear gravity inversion in spherical coordinates with application to the South American Moho

by [Leonardo Uieda](http://www.leouieda.com)
and
[Valéria C. F. Barbosa](http://lattes.cnpq.br/0391036221142471)

**Accepted for publication in the Geophysical Journal International**.

Click on this button to run the code online: [![run on Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/pinga-lab/paper-moho-inversion-tesseroids)

This repository is archived on [figshare](http://figshare.com/): [doi:10.6084/m9.figshare.3987267](https://dx.doi.org/10.6084/m9.figshare.3987267)

A PDF of the article is available at [leouieda.com/papers/paper-moho-inversion-tesseroids-2016.html](http://www.leouieda.com/papers/paper-moho-inversion-tesseroids-2016.html)

Contents:

* Python source code that implements the Moho inversion method proposed in the
  paper (see `code/mohoinv.py`). The modeling code uses the open-source
  geophysics library [Fatiando a Terra](http://www.fatiando.org/).
* Data files required to produce all results and figures presented in the paper.  
* Source code in [Jupyter notebooks](http://jupyter.org/) that produce all
  results and figures (the `.ipynb` files in `code`). You can view static (not
  executable) versions of these notebooks through the links provided in
  `code/README.md`.
* Latex source to produce the manuscript.
* The estimated Moho depth model for South America, available in the file
  `model/south-american-moho.txt` in ASCII xyz format 
  ([download the file](https://raw.githubusercontent.com/pinga-lab/paper-moho-inversion-tesseroids/master/model/south-american-moho.txt)).

The main result from this publication is the gravity-derived Moho depth model
for South America:

![Preview of the estimated Moho depth for South America](https://raw.githubusercontent.com/pinga-lab/paper-moho-inversion-tesseroids/master/model/south-american-moho.jpg?token=AARtIt4v4DyB2aGd81JkbfVlM7sbFqq5ks5W_ClzwA%3D%3D)


## Reproducing the results


### Executing the code online

You can run the Jupyter notebooks online without installing anything
through the free [Binder](http://mybinder.org/) webservice:

[![run on Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/pinga-lab/paper-moho-inversion-tesseroids)

Start by navigating to the `code` folder and choose a notebook to execute.
Beware that the CRUST1.0 and South American Moho notebooks will probably not
be able to run on Binder because they take too long (a few days).
If you want to try out the method start with the `code/synthetic-simple.ipynb`
notebook or see how the figures for the paper are generated in the
`code/paper-figures.ipynb` notebook.
**Tip**: use `Shift + Enter` to execute a code cell.

If you want to run the full inversion on your machine, follow the steps below.


### Running things on your machine

First, download a copy of all files in this repository:

* [Download a zip archive](https://github.com/pinga-lab/paper-moho-inversion-tesseroids/archive/master.zip)
* Or clone the [git](https://git-scm.com/) repository:

        git clone https://github.com/pinga-lab/paper-moho-inversion-tesseroids.git


### Setting up your environment

You'll need a working Python **2.7** environment with all the standard
scientific packages installed (numpy, scipy, matplotlib, etc).  The easiest
(and recommended) way to get this is to download and install the
[Anaconda Python distribution](http://continuum.io/downloads#all).
Make sure you get the **Python 2.7** version.  All the packages that you'll
need are specified in the `environmet.yml` file.  You'll also need to install
the latest development version of the
[Fatiando a Terra](http://www.fatiando.org/) library.

Unzip the contents of this repository (if you've downloaded the zip file) and
`cd` into the root of the repository.  You can use `conda` package manager
(included in Anaconda) to create a virtual environment with all the required
packages installed (including Fatiando). Run the following command in the
repository folder (where `environment.yml` is located):

    conda env create

To activate the conda environment, run

    source activate moho

or, if you're on Windows,

    activate moho

This will enable the environment for your current terminal session.


### Running the inversions and generating figures

The inversion and figure generation are all run inside Jupyter notebooks.  To
execute them, you must first start the notebook server by going into the
repository folder and running (make sure you have the `conda` environment
enabled first):

    jupyter notebook

This will start the server and open your default web browser to the Jupyter
console interface.  In the page, go into the `code` folder and select the
notebook that you wish to view/run.

Here what a notebook looks like running in Firefox:

![Screenshot of the Jupter notebook](https://raw.githubusercontent.com/pinga-lab/paper-moho-inversion-tesseroids/master/screenshot-jupyter-notebook.png?token=AARtIq6LujCeiRLJLIjqQyqAGnV3KS0aks5W_CY1wA%3D%3D)

The notebook is divided cells (some have text while other have code).  Each can
be executed using `Shift + Enter`. Executing text cells does nothing and
executing code cells runs the code and produces it's output.  To execute the
whole notebook, run all cells in order.

All results figures in the paper are generated by the
`code/paper-figures.ipynb` notebook.


## License

All source code is made available under a BSD 3-clause license.  You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors.  See `LICENSE.md` for the full license text.

The model results in file `model/south-american-moho.txt` are available under
the [Creative Commons Attribution 4.0 License (CC-BY)](https://creativecommons.org/licenses/by/4.0/).

The manuscript text is not open source. The authors reserve the rights to the
article content, which has been accepted for publication in the Geophysical
Journal International.
