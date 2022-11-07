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


