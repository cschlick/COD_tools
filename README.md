# COD_tools
A set of scripts to filter the [COD](http://www.crystallography.net/cod/) database for use as a molecular training set.



## Usage

### Get the full dataset using rsync
```console
mkdir COD
cd COD
mkdir -p cif
rsync -av --delete rsync://www.crystallography.net/cif/ cif/
```

This will result in a large number of structure files in the .cif format. These files contain information about the xray experiment, information about the atom coordinates, and information about the bonded atom pairs. 





