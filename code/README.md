# Source code for running the inversions, processing data, and making figures

All data processing, inversion, and figures in this publication are made using
the Jupyter notebooks in this folder.
Follow the links below to see static (not executable) versions of the notebooks.

The Python code that implements the Moho inversion method proposed in the paper
is in file `mohoinv.py`. The inversion code uses the open-source geophysical
modeling library [Fatiando a Terra](http://www.fatiando.org/) version 0.4.
See the following notebooks for examples of using this code:

* [moho-inversion-example.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/moho-inversion-example.ipynb):
  Example of using the Moho inversion class defined in `mohoinv.py`.
* [tesseroid-relief-example.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/tesseroid-relief-example.ipynb):
  Example of using the `TesseroidsRelief` class defined in `mohoinv.py`.

Auxiliary functions and classes for loading and treating the data sets are
defined in `datasets.py`. Examples of using this module:

* [datasets-example.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/datasets-example.ipynb):
  Examples of using the data loading functions in `datasets.py` and the
  class to access the CRUST1.0 model.

The following notebooks produce all results presented in the paper:

* [synthetic-simple.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/synthetic-simple.ipynb):
  Generate the data and run the inversion for the simple synthetic data test.
* [synthetic-crust1.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/synthetic-crust1.ipynb):
  Generate the data and run the inversion for the CRUST1.0 based synthetic data
  test.
* [process-sam-gravity-data.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/process-sam-gravity-data.ipynb):
  Process the raw gravity data of South America: calculate the gravity
  disturbance, terrain correction, and remove the gravitational effects of
  sediments.
* [south-america-moho.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/south-america-moho.ipynb):
  Run the inversion to estimate the Moho depth for South America.

The `results` folder holds the inversion results in zipped Python pickle files.
These files are not meant for user consumption (see the `model` and `data`
folders for more user-friendly formats).
They are loaded in `paper-figures.ipynb` to produce the results figures for the
paper:

* [paper-figures.ipynb](http://nbviewer.jupyter.org/github/pinga-lab/paper-moho-inversion-tesseroids/blob/master/code/paper-figures.ipynb):
  Make all the results figures in the paper.
