# Building TSMP-PDAF Documentation

If you'd like to develop and/or build the TSMP-PDAF documentation, you should:

1. Switch to this directory: `cd docs`
2. Run `pip install -r requirements.txt`
3. (*Optional*) Edit the books source files (`*.md`) in this folder. Check out the [MyST syntax cheat sheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html) for reference.
4. Build the docs: `make clean docs`.
5. (*Optional*) Build TSMP-PDAF source code browser (PDAF-library): `make src-browser-lib`
6. (*Optional*) Build TSMP-PDAF source code browser (TSMP-PDAF interface): `make src-browser-interface`
7. Launch the doc homepage on your default browser: `open _build/html/index.html`

## Contributors

We welcome and recognize all contributions. You can see a list of current contributors in the [contributors tab](https://github.com/HPSCTerrSys/pdaf/graphs/contributors).

## Credits

This project is created using the excellent open source [Jupyter Book project](https://jupyterbook.org/) and the [executablebooks/cookiecutter-jupyter-book template](https://github.com/executablebooks/cookiecutter-jupyter-book).
