for f in *py
do
	python -m pytest --cov-report term-missing --cov=./ $f
done
