# CODE V macro-PLUS programs

This repository provides a collection of CODE V ray-tracing software programs to apply tolerances to the an optical system and retrieve the optical path. 

---

## Getting Started

All functions are located in the `FUNCTIONS` directory. It is necessary to add the path when calling the functions.

Work flow depending on the goal:

- nominal configuration: Nominal_opt.seq and nominalWF.seq;
- optical path: in order to retrieve the optical path run OP_computation;
- perturbations: in order to just apply the perturbations and retrieve the optical path, run in sequence: startup_comp.seq -> data_points49077.seq -> PERTURBATOR.seq;
- compensator: in order to apply the perturbations, apply the compensator and      retrieve the optical path, run in sequence: startup_comp.seq -> data_points49077.seq -> COMPENSATOR.seq; (specify where indicated the compensator to use).

---

## Main Scripts

The main analysis scripts are organized as follows. Each script contains inline comments for further clarification.

### `startup_comp.seq`
- **Description**: initialise functions

### `data_points49077.seq`
- **Description**: set of input data

### `Nominal_opt.seq`
- **Description**: restore nominal optimised configuration

### `nominalWF`
- **Description**: remove surface deformation

### `OP_computation.seq`
- **Description**: compute the optical path

### `normal_distribution.seq`
- **Description**:initialise the variables of the normal distribution

### `tolerwfe.seq`
- **Description**: apply surface deformations

### `PERTURBATOR.seq`
- **Description**: apply tolerances to the system

### `COMPENSATOR.seq`
- **Description**: implement the compensator

### `Control_loop.seq`
- **Description**: compensator based on the bisection method 

### `Opt_comp.seq`
- **Description**: compensator based on the CODE V optimisation commands

---

### External resources

CODE V documentation
