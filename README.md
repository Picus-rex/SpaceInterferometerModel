# Space interferometer modelling

This is a series of functions to study the behaviour of a N apertures nulling interferometer for exoplanets detection.

---

### Description of main files

##### Response functions

This file loads a configuration script of a given array and compute various ideal characteristics, including: 
- The phases optimal splitting following the SVD analysis;
- The classification of the baselines;
- The response function, including, if provided in the configuration files, a simulation of possible disturbances;
- The planet modulation signal over a rotation;
- Several relations about the OPDs and nulling ratios (refer to the functions);
- Perform some plots.

##### Intensity distribution

By loading the default outputs from CODE V, this scripts performs a series of validations of the provided inputs and is currently under definition. 

##### PSF view

This script specifically computes the PSF and some associated properties of a specified array, based on Lay OP. Imaging properties of rotating nulling interferometers (2005). It was mainly used for debug reasons as in the thesis the defined default arrays are considered instead. 

##### Config comparison

This script loads several configurations and computes different ideal performance rates in order to compare different arrays. By default it loads the partial list of Lay's defined configurations with a modified baseline. Eventually, the results are exported in a specified latex table.

##### Generate figures

This scripts exports all the figures accordingly to the specified style for the direct export in the git folder that is synchronised with the thesis repository.


### Description of functions

Please refers to the `docs` folder (currently work in progress). 


---

### External resources

**[yamlmatlab](https://code.google.com/archive/p/yamlmatlab/)** by: 
- Jiri Cigler, Dept. of Control Engineering, CTU Prague http://support.dce.felk.cvut.cz/pub/ciglejir/
- Jan Siroky, Energocentrum PLUS, s.r.o.
- Pavel Tomasko, student at Faculty of Electrical Engineering, CTU Prague and Institute of Chemical Technology, Prague.
MIT License

Adaptation of some functions for external disturbances provided by C. Dandumont.
