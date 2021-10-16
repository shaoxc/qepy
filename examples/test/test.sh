for f in *py
do
	pytest --cov-report term-missing --cov=./ $f
done
