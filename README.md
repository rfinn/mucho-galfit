# mucho-galfit

This repository has the code for running galfit in parallel on WISESize images. 

The galfit models will be used to estimate the size of the star-forming disk from the W3 (12-micron) images. We will compare with the size of the stellar disks, measured from W1 (3.4-micron) images.


# Installation
These instructions are based on these tutorials https://carpentries-incubator.github.io/python-intermediate-development/12-virtual-environments/index.html
## create a virtual environment

```
cd github/mucho-galfit
python3 -m venv venv
```

then activate environment:
```
source venv/bin/activate
```

### If using conda:
```
conda create --name venv
conda activate venv
```

## Install python requirements

```
pip3 install -r requirements.txt
```