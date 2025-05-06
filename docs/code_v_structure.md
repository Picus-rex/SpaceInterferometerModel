# Modelling of N apertures nulling interferometer
### CODE V outputs

Outputs from CODE V must be text files (like .txt) where each simulation is made by $N_P$ rows that follows comma separated values format:

```
OP, X, Y
OP, X, Y
...
```

Each simulation is introduced by a header:

```
== 1
OP, X, Y
OP, X, Y
...
== 2
OP, X, Y
OP, X, Y
...
```

For convenience, the script `add_export_header.py` is provided within the `others` folder to provide a conversion from a file of the shape (no headers, no commas):

```
OP X Y
OP X Y
...
OP X Y
```

### Compensator

In compensators file, the following structure is observed

- Nominal simulation (this is discarded)
- Nominal simulation (with compensator, normally equivalent to nominal)
- Perturbed uncorrected simulation 1
- Perturbed and corrected simulation 1
- Perturbed uncorrected simulation 2
- Perturbed and corrected simulation 2
- Perturbed uncorrected simulation 3
- Perturbed and corrected simulation 3

and so on.
