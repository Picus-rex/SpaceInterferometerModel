# Modelling of N apertures nulling interferometer
### Config file structure

The config file is a yaml file that is parsed by MATLAB using the external library by CTU in Prague and Energocentrum PLUS s.r.o. At the base the file is made by the following structures:

| **Name** | **Description** |
| -------- | --------------- |
| `instrument` | Properties associated to the instrument behaviour, including geometric and leakage properties |
| `environment` | Properties associated to the simulated environment, mainly for leakage computations |
| `simulation` | Properties of the performed simulation |
| `outputs` | Description of outputs to produce |

#### Instrument properties

Fields in italic are optional. This struct contains all the related information on the instrument and some fields may be computed by provided functions, otherwise they can be manually specified.

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `apertures` | -     | Number of apertures of the nulling interferometers |
| `apertures_ratio` | - | Ratio among the maximum and minimum baselines of the system^[In general, it depends on the selected array] | 
| `array` | string    | Name of the used array, either `Diamond`, `Linear` or `X-Array` |
| `baseline` | m      | Maximum baseline of the system |
| *`combination`* | - | Vector of intersity percentage to take from each apertures. Can be derived by other functions. | 
| `efficiencies` | struct | Includes fields `beam_combiner` and `optical_line` as float values |
| `intensities` | -   | Vector of the intensities associated to each aperture |
| *`phase_shifts`* | rad | Can be specified manually or derived using appropriate functions |
| `wavelength` | m    | Wavelength of measurement of the system |

In particular, `combination` and `phase_shifts` are better derived using the provided `compute_optimal_splitting` function but can be manually specified by following the requirements specified on the combinations of apertures conditions to get a nulling interferometer.

#### Environment properties

Fields in italic are optional.

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| *`disturbances`* | struct | See below. All these fields are optional. |
| `exoplanet_position` | rad | Position of the exoplanet in the sky, given as a vector of x and y coordinates |
| `exoplanet_radius` | m | Radius of the exoplanet |
| *`star_radius`* | m   | Radius of the star. Can be computed provided angular radius and distances are provided by other functions. |
| `star_temperature` | K | Star effective temperature. |
| `stellar_angular_radius` | rad | Angular radius of the star |
| `target_distance` | pc | Distance to the study element. |

Where `disturbances` includes:

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `star_flux` | photons/s/m² | Stellar flux received from the star |
| `planet_flux` | photons/s/m² | Flux received from the exoplanet |
| `effective_solid_angle` | sterad | Effective solid angle of the observation |
| *`external_perturbations`* | struct | See below. |
| *`instrumental_leakage`* | struct | See below. |

Where `external_perturbations` includes:

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `intensity` | - | Perturbation intensity |
| `phase` | - | Perturbation phase |
| `x_position` | - | Perturbation x position |
| `y_position` | - | Perturbation y position |

Where `instrumental_leakage` includes:

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `intensity` | -     | 1-Sigma error associated to the aperture |
| `phase`     | -     | 1-Sigma error associated to the phase    |
| `polarisation` | -  | 1-Sigma error associated to the polarisation |

#### Simulation properties

Fields in italic are optional.

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `angular_extension` | multiple | Vector of start, end points [rad] and number of points to consider. Otherwise, is a struct, as described below. |
| `code_v_opd_file`   | -  | See below for fields. |
| `consider_non_ideal` | - | Flag indicating whether to consider non-ideal conditions |
| *`monte_carlo_iterations`* | - | Number of iterations for the Monte Carlo simulation |

If non square arrays are desired (to increase the angular resolution on the x axis for linear arrays, for example, without exceeding resources), then `angular_extension` must be an array with the following fields:

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `theta_x`| multiple | Vector of start, end points [rad] and number of points to consider for the x axis. |
| `theta_y`| multiple | Vector of start, end points [rad] and number of points to consider for the y axis. |

Where `code_v_opd_file` can include the following fields; if specified, the additional field `code_v` will be present once `convert_data` is runned (notice it can slow down computations): 

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `compensator` | -   | Path to compensator export from CODE V (*). |
| `nominal`     | -   | Path to nominal export from CODE V.   |
| `perturbed` | -     | Path to perturbed export from CODE V. |

(*) If compensator is provided, all the other fields are always overwritten. For the structure of CODE V export, see the [respective file](code_v_structure.md). 

#### Output properties

Fields in italic are optional.

| **Name** | **Unit** | **Description** |
| -------- | -------- | --------------- |
| `plot_array` | - | Flag indicating whether to plot the array |
| `plot_planets` | - | Flag indicating whether to plot the planets (0 means no plot, 1 means simplified planes, 2 means random exoplanet from P-Pop) |
| `planets_angular_separations` | m | List of angular separations for the planets (used for `plot_planets` set to 1). |
| `plot_star` | - | Flag indicating whether to plot the star |
