setup:
	conda create -n mohoinv --yes pip python2.7
	source activate mohoinv
	conda install --yes --file requirements.txt
	pip install https://github.com/fatiando/fatiando/archive/inversion.zip

clean:
	rm -rf notebooks/*.png notebooks/profiling.txt
	find . -name "*.pyc" -exec rm -v {} \;
