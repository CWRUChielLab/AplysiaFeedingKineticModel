# *Aplysia* Feeding Kinetic Model

Code repository: https://github.com/CWRUChielLab/AplysiaFeedingKineticModel

This repository contains code for a detailed biomechanical model of the feeding
apparatus of the sea slug *Aplysia californica*. See Snyder (2005) for a detailed
"Guidebook to the Model" defining equations and parameters:

> Snyder, V. A. (2005). Analysis of the biomechanics and neural control of two kinetic models of the buccal mass of Aplysia [M.S. thesis, Case Western Reserve University]. http://search.proquest.com/docview/305390183/abstract/B303E6FBC4924763PQ/1

The model is instantiated in C++, and visualization scripts are written in
Python.

## Usage

1. Navigate to the project directory:
```
cd AplysiaFeedingKineticModel
```

2. The [`Makefile`](./Makefile) provides convenient recipes for compiling the
code, batch-running the model, and validating results. To simply compile the
source code, run the following:
```
make
```
Optionally, use `make debug` to compile the code with extra print commands for
debugging. The `make clean` command can be used to clean up your working
directory by removing compile artifacts and model outputs.

3. Run the model executable:
```
./model <behavior>
```
where `<behavior>` is one of several valid behavior types that can be
simulated, such as `Bite`. For a full list of possible behavior types, run
`./model` without arguments. Each type differs in the activation timing of
muscles, resulting in different movements.

Available behavior types:
* `Bite`: A bite behavior (failure to grasp food)
* `SwallowA`: A type-A swallow (smaller amplitude)
* `SwallowB`: A type-B swallow (larger amplitude with dorsal rotation of the held food)
* `SwallowPerturbed`: A modified type-B swallow with a brief external force applied to the food
* `RejectionA`: A type-A rejection (smaller amplitude)
* `RejectionB`: A type-B rejection (larger amplitude with ventral rotation of the held object)
* `IzhikevichBite`: Similar to `Bite`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzhikevichSwallowA`: Similar to `SwallowA`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzhikevichSwallowB`: Similar to `SwallowB`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzhikevichSwallowPerturbed`: Similar to `SwallowPerturbed`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzhikevichRejectionA`: Similar to `RejectionA`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzhikevichRejectionB`: Similar to `RejectionB`, but with a layer of integrate-and-fire motor neurons driving the muscles at defined frequencies
* `IzExampleSwallow`: A more advanced model of the neuronal inputs to the muscles, where many identified motor neurons are represented
* `IzSwallowBmoreaccurate`: Similar to `IzExampleSwallow` (?)
* `Dynamic`: A version that uses a user-customizable input file (CSV format) to specify the timing and magnitude of inputs to the muscles
* `NeuromechanicalInput`: Similar to `Dynamic`, but expects a different set of inputs created by another model of the motor neurons with interneurons coordinating them (https://github.com/CWRUChielLab/AplysiaNet)

Integrate-and-fire neurons were implemented following Izhikevich (2003):

> Izhikevich, E. M. (2003). Simple model of spiking neurons. IEEE Transactions on Neural Networks, 14(6), 1569–1572. https://doi.org/10.1109/TNN.2003.820440

## Visualizing Results

When the model runs, it generates output files in CSV format containing
periodic snapshots of the values of various quantities, including positions,
angles, and levels of muscle activation. Python scripts can generate plots and
animations from these output files:
* `newanimation.py`: generates movies of the motions of the buccal mass, with time plots of the muscle inputs
* `PlotVariablesSlugOutput.py`: generates time plots of the kinematic variables
* `makerasterplot.py` and `PlotVariablesIzhikevich.py`: generate figures of the activity of the integrate-and-fire motor neurons, which are not included in all behavior types

The `make animations` and `make figures` commands can be used to batch-run the
model on a subset of the behavior types, renaming and saving the resulting
movies and figures to dedicated directories.

The Python scripts require a handful of packages that can be installed with
`pip`:

    pip install numpy pandas matplotlib tqdm

## Validation

The outputs of the model when run on a subset of the behavior types are saved
in the [`Check`](./Check) directory. The `make check` command can be used to
systematically compare your version's outputs to these "gold standards". This
is useful for catching differences in implementation of libraries across
platforms, or other discrepancies that might cause the model to run differently
across systems. Deviations should be investigated as they may lead to bugs.

Of course, if the code is modified intentionally in such a way that the outputs
change (e.g., a new column is added to the output files), the "gold standard"
against which `make check` runs will need to be updated. The `make bless`
command will systematically run the model on the subset of behavior types and
copy the outputs to the [`Check`](./Check) directory. These changed outputs,
once verified to be correct, should be committed to the repository so that
results can be validated against them in the future (they will become the new
"blessed" gold standard of what outputs the model _should_ produce).
