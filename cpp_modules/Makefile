export PYTHON_VERSION := 3.12

# sudo apt install python$(PYTHON_VERSION)-full python$(PYTHON_VERSION)-dev libgudhi-dev

# deactivate
# python$(PYTHON_VERSION) -m venv ~/venv/tda
# . ~/venv/tda/bin/activate
# pip install --upgrade pip
# pip install pipreqs
# pipreqs --force <project root path>
# pip install -r requirements.txt
# python -m pip install pybind11

BLDDIR = build

all release debug test:
	$(MAKE) -C adrcm $@
	$(MAKE) -C hypergraph $@
	$(MAKE) -C network $@

clean :
	@rm -rf *~ core $(BLDDIR)