# kappamin
A python3 code for calculations of the minimum limit to thermal conductivity

## Features

- Models of the minimum limit to thermal conductivity under Cahill assumption[^1]
  - Debye model[^2]
  - BvK (Born–von Karman) model[^3]
  - Pei model[^4]
- Temperature-dependence
  - Finite temperature
  - Ideal infinite temperature
- Relative
  - Heat Capacity
  - Minimum mean-free-path
  - Minimum average phonon lifetime
- Running mode
  - Command line mode based on a simple configuration file (for the routine analysis)
  - Based on prepared scripts (for general researchers without programming skills)
  - Use as a python module (for expert usage)

## Getting Started

In order to install the module, you need to install the dependencies first:

- python3
- numpy
- scipy

Download the package from GitHub website or using `git`,

```
$ git clone https://github.com/JianboHIT/kappamin.git
```

then run `setup.py` to install.

```
$ cd kappamin
$ python3 setup.py install
```

Python module `kappamin` will be installed if no error. 
In order to invoke `kappamin` module, you need to prepare a configuration file 
(see [Example_Debye.txt](Example_Debye.txt)
and [Example_BvK.txt](Example_Debye.txt)
in the source package).

```
$ python -m kappamin [KAPPAMIN.txt]
```

Here `KAPPAMIN.txt` indicates the filename of configuration file. 
It is worth mentioning that the filename is optional. 
If the filename is not given, the program will read the file named as `KAPPAMIN.txt` if it existed.

Alternately, a more convenient way to implement calculation is by an executable script
(see [ExceuteScript.py](ExecuteScript.py)), then run it by python3. 
On Linux or Windows Terminal:

```
$ python ExceuteScript.py
```

On Windows, if it has been configured that the default program to open .py file is python3,
you just need to move ExceuteScript.py to the directory at where the configuration file is located
and double-click it to run.

## Feedback and report bugs

See [GitHub Issue page](../../issues).

## Change log

(More details see [CHANGELOG](CHANGELOG))

- 2022.10.06 v0.0.1 Initial package version
- 2022.10.16 v0.1.0 Develop Debye, BvK, and Pei models


<br/><br/>

## Reference

[^1]: D.G. Cahill, R.O. Pohl, Heat flow and lattice vibrations in glasses, Solid State Communications, 70 (10) (1989) 927-930. [https://doi.org/10.1016/0038-1098(89)90630-3](https://doi.org/10.1016/0038-1098(89)90630-3)

[^2]: P. Debye, Zur theorie der spezifischen wärmen, Annalen Der Physik, 344 (14) (1912) 789-839. [https://doi.org/10.1002/andp.19123441404](https://doi.org/10.1002/andp.19123441404)

[^3]: M. Born, T. Von Karman, Vibrations in space gratings (molecular frequencies), Z Phys, 13 (1912) 297-309.

[^4]: Z. Chen, X. Zhang, S. Lin, L. Chen, Y. Pei, Rationalizing phonon dispersion for lattice thermal conductivity of solids, National Science Review, 5 (6) (2018) 888-894. [https://doi.org/10.1093/nsr/nwy097](https://doi.org/10.1093/nsr/nwy097)
