# Space interferometer modelling

This repository contains a series of functions to study the behaviour of a N apertures nulling interferometer for exoplanets detection from a numerical point of view, looking at important metrics (nulling ratio, modulation efficiency, exoplanets yield), including the study of optical perturbations by post-processing outputs from ray-tracing softwares.

---

### Getting Started

The function library is provided within the `function` folder. They can be called individually or exploited like shown in the main files within the root of the folder. All the main files starts by clearing the workspace and by loading the functions folders, therefore everything works out of the box, provided that the necessary toolboxes are installed.

Most of the scripts (whose description is provided below) rely on configuration files, written in YAML. Examples are given within the `config` folder. The documentation folder `docs` contains a document that explains the main elements of a configuration file. Some of them (the ones that concern perturbations), also require external outputs from ray-tracing software (more about them in the [CODE V integration](docs/code_v_structure.md)): outputs from that files are given, as an example, in the `code_v` folder. 


### Description of main files

The complexity associated to each task is here presented in a sequential way, therefore the former files execute simpler tasks with respect to the latters. Each file might include a more detailed explanation for each action. 

##### Response functions
*Requires:* configuration file.

This file loads a configuration script of a given array and compute various ideal characteristics, including: 
- The phases optimal splitting following the SVD analysis;
- The classification of the baselines;
- The response function, including, if provided in the configuration files, a simulation of possible disturbances;
- The planet modulation signal over a rotation;
- Several relations about the OPDs and nulling ratios (refer to the functions);
- Perform some plots.

##### PSF view
*Requires:* configuration file.

This script specifically computes the PSF and some associated properties of a specified array, based on Lay OP. Imaging properties of rotating nulling interferometers (2005). It was mainly used for debug reasons as in the thesis the defined default arrays are considered instead. 

##### Config comparison
*Requires:* configuration file.

This script loads several configurations and computes different ideal performance rates in order to compare different arrays. By default it loads the partial list of Lay's defined configurations with a modified baseline or the most available data on LIFE, Darwin and TPF-I. Eventually, the results are exported in a specified latex table.

##### Tolerancing analysis
*Requires:* configuration file and optical data.

This script analyses the optical data in an aggregate way to compute some figures that study how perturbation affects the system.

##### Exoplanets yields
*Requires:* configuration file and optical data.

This script computes the exoplanet yield both for the nominal system and the perturbed system (using aggregated data). It can be partially executed with the configuration file alone by stopping early in the computation.

##### Config comparison perturbated
*Requires:* configuration file and optical data.

It repeats what it has been done in `config_comparison` but including some perturbations (it takes longer to run), presenting some aggregated results.

##### Interferogram arms
*Requires:* configuration file and optical data.

It computes some perturbation effects at the pupil plane, therefore characterising the effects of the interferometer perturbations for every propagated points in the analysis. 

##### Generate figures
*Requires:* configuration file and optical data.

This scripts exports all the figures accordingly to the specified style for the direct export in the git folder that is synchronised with the thesis repository. This script may not be of interest in the global analysis and can be skipped. 

### Description of functions

Each function (at least those that have a larger impact on the code base) has its own `help`. By running `help function_name`, simple usage and input/output explanation is provided. 

---

### External resources

This project includes, with their respective licenses, the following elements:

**[yamlmatlab](https://code.google.com/archive/p/yamlmatlab/)** by: 
- Jiri Cigler, Dept. of Control Engineering, CTU Prague http://support.dce.felk.cvut.cz/pub/ciglejir/
- Jan Siroky, Energocentrum PLUS, s.r.o.
- Pavel Tomasko, student at Faculty of Electrical Engineering, CTU Prague and Institute of Chemical Technology, Prague.
MIT License

Adaptation of some functions for external disturbances provided by C. Dandumont.
