# G4Hunter: Transfer to Python3

([Re-evaluation of G-quadruplex propensity with G4Hunter](http://nar.oxfordjournals.org/content/early/2016/01/19/nar.gkw006.full.pdf+html))

- Transfer the origin script in python2 version to python3 version
- Remove the visualization function(If you need it, you can modify the script referencing the py2 script)

## Requirements & Obtaining

- **Python3.7** is recommended

- **Biopython**: version 1.76

- **Numpy**: version 1.16.5

  ```shell
  # To install Bio, numpy:
  pip install biopython
  pip install numpy
  ```



## Launching

```shell
python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>
```

