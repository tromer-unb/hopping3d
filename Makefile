PYTHON ?= python3
.PHONY: install test examples figures reproduce-full manuscript clean
install:
	$(PYTHON) -m pip install -e .
test:
	pytest -q
examples:
	hopping3d examples/first_passage/params.json
	hopping3d examples/diffusion_tensor/params.json
figures:
	$(PYTHON) reproduction/scripts/figure1_overview.py
	$(PYTHON) reproduction/scripts/plot_publication_figures.py
	$(PYTHON) reproduction/scripts/figureS_controls.py
	$(PYTHON) reproduction/scripts/figureS_convergence.py
reproduce-full:
	bash reproduction/reproduce.sh full
manuscript: figures
	bash manuscript/build.sh
clean:
	rm -rf reproduction/generated/* examples/*/results .pytest_cache build dist *.egg-info
