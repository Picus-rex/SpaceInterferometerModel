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

- Perturbed uncorrected simulation 1
- Perturbed and corrected simulation 1
- Perturbed uncorrected simulation 2
- Perturbed and corrected simulation 2
- Perturbed uncorrected simulation 3
- Perturbed and corrected simulation 3

and so on, therefore, when the `compensator` field is provided, the field `perturbed` is overwritten (a warning will be displayed). 

### Splitted files 

Since GitHub does not allow for larger files to be upload, it is possible to load splitted files. Splitted files, generated with `split_large_files.py` function, follow the form `basename.id.txt`, where `id` starts from 1 up to all the necessary files. So that, in the configuration file, one can simply indicate, for example:

```
compensator: "code_v/0510_compensator_bis.txt"
```

Then all the files `code_v/0510_compensator_bis.1.txt`, `code_v/0510_compensator_bis.2.txt` etc. will be imported. 
