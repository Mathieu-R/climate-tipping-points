# climate-tipping-points

Investigating tipping points in climatology. In particular, a fold-hopf bifurcation.    
This work has been based on this article from Dekker et al. : https://esd.copernicus.org/articles/9/1243/2018/  

![bifurcations](bifurcations.png)

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

Install required packages
```bash
$ python3 -m pip install numpy matplotlib scipy click tqdm
```

Launch examples
```bash
$ python3 index.py --plot=time-series 
$ python3 index.py --plot=bifurcations
```
