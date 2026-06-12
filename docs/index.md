# airsspy Documentation

Modern Python interface for AIRSS (Ab initio Random Structure Searching).

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Getting Started
:link: getting-started/index
:link-type: doc

New to airsspy? Start here to learn installation and basic usage.
:::

:::{grid-item-card} How-To Guides
:link: how-to/index
:link-type: doc

Step-by-step guides for common tasks and workflows.
:::

:::{grid-item-card} Explanations
:link: explanation/index
:link-type: doc

Understand AIRSS concepts, seed files, and buildcell parameters.
:::

:::{grid-item-card} API Reference
:link: reference/index
:link-type: doc

Complete API documentation for all modules and classes.
:::

:::{grid-item-card} Command Line
:link: reference/cli
:link-type: doc

Reference for the installed `ap` command and local AIRSS workflows.
:::
::::

## Quick Example

```python
from airsspy import SeedAtoms

# Create a seed for silicon structure search
seed = SeedAtoms('Si8', cell=[5, 5, 5], pbc=True)
seed.gentags.minsep = 2.0
seed.gentags.varvol = 20

# Generate a random structure
random_atoms = seed.build_random_atoms()
```

## What is airsspy?

airsspy provides a modern, Pythonic interface to AIRSS by:
- Extending ASE's Atoms class with AIRSS-specific features
- Managing buildcell subprocess communication
- Handling AIRSS/CASTEP `.res` file I/O and lossless RES/extxyz conversion
- Running local searches, relaxations, single-points, and packed-file ranking
- Deploying jobflow searches and querying stored results

[Learn more →](explanation/airss-overview.md)

## Features

- **ASE Integration**: Seamless integration with the Atomic Simulation Environment
- **Modern Python**: Type hints, clean API, Python 3.9+ support
- **Structure Generation**: Generate random structures with flexible constraints
- **RES File Support**: Read and write CASTEP .res files
- **Buildcell Interface**: Full access to AIRSS buildcell parameters
- **CLI Workflows**: `ap run`, `ap rank`, `ap convert`, `ap pack`, and jobflow database tools

## Installation

Install airsspy using pip:

```bash
pip install airsspy
```

For development installation with all dependencies:

```bash
git clone https://github.com/zhubonan/airsspy.git
cd airsspy
pip install -e ".[dev,docs]"
```

See the [Installation Guide](getting-started/installation.md) for more details.

## Quick Links

- [Quickstart Tutorial](getting-started/quickstart.md) - Your first structure search
- [Creating Seeds](how-to/create-seeds.md) - How to create structure templates
- [Command Line](reference/cli.md) - Installed `ap` commands and examples
- [API Reference](reference/index.md) - Complete API documentation
- [GitHub Repository](https://github.com/zhubonan/airsspy) - Source code and issues

```{toctree}
:hidden:

getting-started/index
how-to/index
explanation/index
reference/index
```
