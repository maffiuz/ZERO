# ZERO

> A zero-dimensional, lumped parameter model, capable of solving a simplified nuclear fusion plasma energy balance equation.

![Python](https://img.shields.io/badge/python-3.14%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

This tool was initially developed as a Bachelor Thesis by the Author at [Politecnico di Torino](https://www.polito.it/). Then it was expanded, to be a suitable demonstrative teaching tool, capable of showing the resolution of a lumped-parameter energy balance equation of the plasma in a nuclear fusion reactor.

## Features

- Numerical solution of the time dependent energy balance equation of a nuclear fusion plasma
- Numerical solution of the evolution of the plasma composition, in form of particle densities of deuterium, tritium, ashes, and impurities
- Possibility of selecting specific nuclear reactions, operative conditions and impurity puffing
- The code produces plots that are explained and can be easily consulted thanks to the built-in Jupyter notebook

## Installation

Requires Python 3.14+. Download the repository or use:

```bash
git clone https://github.com/maffiuz/ZERO.git
cd ZERO/ReactorSimulator
```

The virtual environment is based on the *pipenv* library:

```bash
pip install pipenv
pipenv shell
```

Build the dependency refences, based on the *Pipfile* with:
```
pipenv lock
```

Download all the dependencies with:

```
pipenv sync
```

## Atomic data

The code is based on atomic data from the Atomic Data and Structure Analysis (ADAS) database, particularly through the OPEN-ADAS website: https://open.adas.ac.uk . Currently the data has to be downloaded manually. For further information, contact the Author.

## Usage

### CLI
First activate the virtual environment, then run the jupyter notebook:

```bash
pipenv shell
jupyter notebook simulator_report.ipynb​
```

## Configuration

The configuration can be done in two ways:
1. By using a **run_report_** file, copyied into the *ReactorSimulator* folder
2. By modifying variables directly in the source code, in **userDefinded/user_inputs.py**

In presence of a report, the program will always ask the user whether to use it or the data in the source.

Here is a list of the main input parameters:

| Parameter | Name in *run_report* | Description | Units |
|-----------|----------------------|-------------|-------|
| `t_start` | *t_start* | simulation starting time | [s] |
| `t_end` | *t_end* | simulation ending time | [s] |
| `n_points` | *n_points* | number of points on which the solver projects the solution | [-] |
| `a` | *a* | plasma minor radius | [m] |
| `R` | *R* | plasma major radius | [m] |
| `t_ref` | *t_ref* | refueling starting time | [s] |
| `D0_in` | *D0_in* | deuterium refueling rate | [s-1] |
| `T0_in` | *T0_in* | tritium refueling rate | [s-1] |
| `TimeDown` | *t_down* | startup time | [s] |
| `Aux_heat_start` | *W_ext* | auxiliary heating during startup | [W] |
| `Aux_heat_operat` | *W_ext* | auxiliary heating | [W] |
| `IP` | *IP* | plasma current | [A] |
| `BT` | *BT* | toroidal magnetic field component | [T] |
| `E_START` | *E* | total initial energy content | [J] |
| `ND0_START` | *nD0* |  initial neutral deuterium particle density | [m-3] |
| `ND1_START` | *nD1* | initial ionized deuterium particle density | [m-3] |
| `NT0_START` | *nT0* | initial neutral tritium particle density | [m-3] |
| `NT1_START` | *nT1* | initial ionized tritium particle density | [m-3] |
| `NHE4_START` | *nHe4_0* | initial neutral He4 particle density | [m-3] |
| `NHE41_START` | *nHe4_1* | initial ionized He4 particle density | [m-3] |
| `NHE42_START` | *nHe4_2* | initial alpha particle density | [m-3] |
| `IMPURITIES` | Multiple | starting density, injection rate, injection start time, injection time | |
| `REACTIONS` | Multiple | boolean switches to activate/deactivate nuclear fusion reactions | |


## Roadmap

- [ ] Interface to download ADAS data


## Development

This code is developed and maintained by R.V. Maffiodo in collaboration with F.Subba (POLITO).

Contributions are welcome. Please open an issue first to discuss changes or contact the Author.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.