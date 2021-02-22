# climate-tipping-points

Investigating tipping points in climatology. In particular, a fold-hopf bifurcation.

### Running simulations
Clone the project
```bash
$ git clone https://github.com/Mathieu-R/climate-tipping-points.git
```

Create virtual environment
```bash
$ python3 -m venv <env-name>
$ source env/bin/activate
$ python3 -m pip install --upgrade pip
```

Connect Jupyter to virtual environment
```bash
$ python3 -m pip3 install ipykernel
$ python3 -m ipykernel install --name=<env-name>
```

Install required packages
```bash
$ python3 -m pip install --requirement=requirements.tex
```

Launch Jupyter Lab
```bash
$ jupyter lab
```

### Using Kite extension with Jupyter Lab
```bash
$ python3 -m pip install "jupyterlab-kite>=2.0.2"
```
