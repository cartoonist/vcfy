PIP=pip
TEST=nosetests
PKG_NAME=vcfy
TESTREPO=pypitest
MAINREPO=pypi

# Specifying phony targets
.PHONY: init test dist-test dist

init:
	${PIP} install -r requirements.txt

test:
	${TEST}

README.rst: README.md
	pandoc -o $@ $<

dist-test: README.rst
	python setup.py register -r ${TESTREPO}
	python setup.py sdist upload -r ${TESTREPO}

dist: README.rst
	python setup.py register -r ${MAINREPO}
	python setup.py sdist upload -r ${MAINREPO}
