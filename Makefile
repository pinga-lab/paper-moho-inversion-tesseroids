clean:
	rm -rf notebooks/*.png notebooks/profiling.txt
	find . -name "*.pyc" -exec rm -v {} \;
