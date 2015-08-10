setup: mkenv install_requires install_fatiando

mkenv:
	conda create -n mohoinv --yes pip python=2.7

install_requires:
	bash -c "source activate mohoinv && conda install --yes --file requirements.txt"

install_fatiando:
	bash -c "source activate mohoinv && pip install --upgrade https://github.com/fatiando/fatiando/archive/inversion.zip"

delete_env:
	bash -c "source deactivate; conda env remove --name mohoinv"

clean:
	rm -rf notebooks/*.png notebooks/profiling.txt
	find . -name "*.pyc" -exec rm -v {} \;
