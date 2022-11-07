# COD Tools
A set of scripts to filter the [COD](http://www.crystallography.net/cod/) database for use as a molecular training set.


## Installation
Install additional dependencies into the same environment
```console
conda install -c conda-forge tqdm rdkit networkx cctbx-base mdanalysis
```

## Usage

### Get the full dataset using rsync
```console
mkdir COD
cd COD
mkdir -p cif
rsync -av --delete rsync://www.crystallography.net/cif/ cif/
```

This will take a while. It downloads a large number of small molecule structure files in a .cif text format. These files contain information about the xray experiment and data, information about the atom coordinates, and information about the bonded atom pairs. 





