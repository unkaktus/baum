
### BAM integration tests

First, install MambaForge (https://mamba.readthedocs.io/en/latest/mamba-installation.html).

All the actions below are performed from `ci` directory.

Then, create an environment named `bam-ci` from `environment.yaml` file

```shell
mamba env create -n bam-ci -f environment.yaml
```

This will install all the required packages for tests, including `task` command.

Activate the `bam-ci` environment:

```shell
mamba activate bam-ci
```

Build and install XDMF library package:

```shell
task install-xdmf
```

For the test with tabulated EOS, you need to generate the tables
from CompOSE data, first installing `compose`:

```shell
task install-compose
task generate_EOS_table_SFHo
```

Build BAM:

```shell
task build
```

Run all the tests:

```shell
task run-all
```

# Optional

If you wish to run different tests separately, list all available tasks via:

```shell
task --list-all
```