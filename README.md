# Space Interferometer Modelling

This repository provides a collection of MATLAB functions for numerically studying the behavior of a nulling interferometer with *N* apertures for exoplanet detection. It evaluates key performance metrics such as nulling ratio, modulation efficiency, and exoplanet yield. The library also supports the analysis of optical perturbations through post-processing of outputs from ray-tracing software (e.g., CODE V).

The repository is modular and intended for both performance analysis and simulation-based tolerancing studies. It is structured to support a variety of use cases, from ideal response modelling to perturbed system analysis.

---

## Getting Started

All functions are located in the `functions` directory. They can be called individually or used as shown in the example scripts in the root folder. These scripts automatically clear the workspace and add the necessary paths, so they work out of the box, provided that the required MATLAB toolboxes are installed.

> ⚠️ Some scripts expect specific folders to exist. If a `File not found` error occurs, create the missing folders or update the script paths accordingly. In particular, the `export` folder must exist in the root directory.

Most scripts rely on configuration files written in YAML format. Examples are available in the `config` folder. Documentation describing the structure and syntax of the configuration files can be found in the `docs` folder. Some analyses (e.g., those involving perturbations) also require ray-tracing outputs, which are expected to be structured as described in [docs/code_v_structure.md](docs/code_v_structure.md). Sample outputs are provided in the `code_v` folder.

> ⚠️ Some signal modelling functions (e.g., `add_external_sensitivity`) depend on proprietary software and are not included in this repository.

---

## Main Scripts

The main analysis scripts are organized in increasing complexity. Each script contains inline comments for further clarification.

### `response_functions.m`
- **Requires**: Configuration file
- **Description**: 
  - Loads a configuration for a specific interferometer array
  - Computes ideal system characteristics including:
    - Optimal phase splitting via SVD
    - Baseline classification
    - Response function (optionally with disturbances)
    - Planet modulation signal over rotation
    - Relations between OPDs and nulling ratios
  - Generates several plots for visualization

### `psf_view.m`
- **Requires**: Configuration file
- **Description**: Computes the point spread function (PSF) and related imaging properties based on Lay OP. Imaging properties of rotating nulling interferometers (2005). Mainly used for debugging.

### `config_comparison.m`
- **Requires**: Configuration file
- **Description**: Compares the ideal performance of multiple configurations (e.g., Lay configurations, LIFE, Darwin, TPF-I). Outputs a LaTeX table with the comparison results.

### `tolerancing_analysis.m`
- **Requires**: Configuration file and optical data
- **Description**: Aggregates and analyses optical perturbation data to evaluate their impact on system performance.

### `exoplanets_yields.m`
- **Requires**: Configuration file and optical data (optional for partial execution)
- **Description**: Computes the exoplanet yield for both nominal and perturbed configurations using aggregate performance data.

### `config_comparison_perturbed.m`
- **Requires**: Configuration file and optical data
- **Description**: Similar to `config_comparison`, but includes the effects of perturbations. More computationally intensive.

### `interferogram_arms.m`
- **Requires**: Configuration file and optical data
- **Description**: Propagates perturbations at the pupil level to evaluate their spatial effects across the interferometer arms.

### `generate_figures.m`
- **Requires**: Configuration file and optical data
- **Description**: Exports all relevant figures in a predefined format for use in a thesis repository. Optional for general users.

---

## Function Documentation

All major functions include inline documentation. Run the following command in MATLAB for usage information:

```matlab
help function_name
```

This includes descriptions of inputs, outputs, and example usage.

---

### External resources

This project includes, with their respective licenses, the following elements:

**[yamlmatlab](https://code.google.com/archive/p/yamlmatlab/)** by: 
- Jiri Cigler, Dept. of Control Engineering, CTU Prague http://support.dce.felk.cvut.cz/pub/ciglejir/
- Jan Siroky, Energocentrum PLUS, s.r.o.
- Pavel Tomasko, student at Faculty of Electrical Engineering, CTU Prague and Institute of Chemical Technology, Prague.
MIT License
