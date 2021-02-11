# Get Started

## Content

- [Run dfttk with `dfttk run` command](#Run-dfttk-with-dfttk-run-command)
  - [Run dfttk for a single structure or a single folder](#Run-dfttk-for-a-single-structure-or-a-single-folder)
  - [Write out the workflow](#Write-out-the-workflow)
  - [Run dfttk for a folder and its sub-folder](#Run-dfttk-for-a-folder-and-its-sub-folder)
  - [Run dfttk with individual settings](#Run-dfttk-with-individual-settings)
  - [Submit job to launchpad and queue system](#Submit-job-to-launchpad-and-queue-system)
- [Run dfttk with python script](#Run-dfttk-with-python-script)
- [Help for `dfttk run` command](#Help-for-dfttk-run-command)

---

## Run dfttk with `dfttk run` command

```shell
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-l] [-m [MAX_JOB]] [-o]
```

The help of `dfttk run` command, please ref. [Help for `dfttk run` command](#Help-for-dfttk-run-command) or simply run `dfttk run -h`

All the examples is run under current folder (get_started)

```shell
get_started
 NaCl.cif
 NotSupportFile
 POSCAR
 test_folder
     Fe3Ni.cif
     Fe3Ni-settings.yaml
     POSCAR
     POSCAR-2
     SETTINGS-POSCAR.yaml
     SETTINGS.yaml
     sub_folder
         POSCAR
```

### Run dfttk for a single structure or a single folder

`-f` parameter, specify the structure or folder containing structure

```shell
dfttk run -f STRUCTURE/STRUCTURE_FOLDER
```

- When the parameter value of `-f` is a file  it run this file; if it is a folder, it will search all supported files in the specified folder.
- It support **POSCAR**, **CONTCAR**, **CHGCAR**, **LOCPOT**, **cif**, **vasprun.xml**, **pymatgens structures**, **mcsqs's structure**, more details, ref [pymatgen.core.structure.IStructure.from_file](https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure.from_file)
  - For **POSCAR** or **CONTCAR**, the name should be **\*POSCAR\*** or **\*CONTCAR\*** or **\*.vasp**
  - For **CHGCAR** or **LOCPOT**, the name should be **\*CHGCAR\*** or **\*LOCPOT\***
  - For **cif**, the name should be **\*.cif*** or **\*.mcif\***
  - For **vasprun.xml**, the name should be **vasprun\*.xml\***
  - For **pymatgen's structure**, the name shold be **\*.json** or **\*.yaml**
  - For **atat's structure**, the name should be **\*rndstr.in\*** or **\*lat.in\*** or **\*bestsqs\***

[TO TOP](#Content)

### Write out the workflow

`-o` parameter, write out work flow

- When add `-o` parameter, the work flow for every structure will be write out

- The file name of the work flow is `dfttk_wf-` + filename of the structure + `.yaml` (e.g. `dfttk_wf-POSCAR.yaml`)

  ```shell
  dfttk run -f POSCAR -o
  ```

  It will write out the work flow in current folder as `dfttk_wf-POSCAR.yaml`

- The file will write out in the same folder with corresponding structure.

  ```shell
  dfttk run -f ./test_folder/POSCAR -o
  ```

  It will write out the work flow as `./test_folder/dfttk_wf-POSCAR.yaml`

[TO TOP](#Content)

### Run dfttk for a folder and its sub-folder

`-r` parameter, recursive the folder

`-mh` parameter, specify the match pattern

- Add `-r` parameter, it will searching the folder (specified by `-f` parameter) and its sub-folder

  ```shell
  dfttk run -f . -r
  ```

  It will run the following structure

  ```shell
  get_started/NaCl.cif
  get_started/POSCAR
  get_started/test_folder/Fe3Ni.cif
  get_started/test_folder/POSCAR
  get_started/test_folder/POSCAR-2
  get_started/test_folder/sub_folder/POSCAR
  ```

-  Add `-mh` parameter, it will filter the filename by the value of `-mh`, the pattern should be placed in quotes ('' or "")

  ```shell
  dfttk run -r -mh '*.cif'
  ```

  It will run the following structure

  ```shell
  get_started/NaCl.cif
  get_started/test_folder/Fe3Ni.cif
  ```

[TO TOP](#Content)

### Run dfttk with individual settings

`-s` parameter, specify the name for settings, default: SETTINGS

- `SETTINGS.yaml` or `SETTINGS.json`: The global settings for current folder

- `SETTINGS-*.yaml(json)` or `*-SETTINGS.yaml(json)`: The individual settings, `*` is the name of the structure (without ext)

- Case insensitive. (both SETTINGS and settings are OK)

- The value of `-s` parameter will replace SETTINGS.

- Examples:

  ```shell
  dfttk run -f test_folder -s set
  ```

  It will take `./test_folder/Fe3Ni-SET.yaml` as setting file, and `./test_folder/Fe3Ni-settings.yaml` will not taken as setting file.

[TO TOP](#Content)

### Submit job to launchpad and queue system

`-l` parameter, launch to launchpad

`-m` parameter, launch to queue system and determine the number of jobs

- Add `-l` parameter, it will submit the job to launchpad

- Add `-m N` parameter, it will submit the job to queue system. `N` specifies the number of jobs run at the same time. (Note: This only work when `-l` parameter is added.) 

  ```shell
  dfttk run -l
  dfttk run -l -m 1
  dfttk run -l -m 2
  ```

  - When `-m 1`, it will run `qlaunch singleshot` (fireworks command)
  - When `-m N` (N>1), it will run `qlaunch rapidfire -m N` (fireworks command)
  - For more details, please ref. [Fireworks: Launch Rockets through a queue](https://materialsproject.github.io/fireworks/queue_tutorial.html)



[TO TOP](#Content)

## Run dfttk with python script





[TO TOP](#Content)

---

## Help for `dfttk run` command

```shell
dfttk run -h
```

```shell
DFTTK version: 0.1+98.ge7aa39c.dirty
Copyright  Phases Research Lab (https://www.phaseslab.com/)

usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -h, --help            show this help message and exit
  -f STRUCTURE_FOLDER, --structure_folder STRUCTURE_FOLDER
                        The folder/file containing the structure, Default: '.'
  -mh MATCH_PATTERN, --match_pattern MATCH_PATTERN
                        The match pattern for structure file, e.g. *POSCAR*,
                        Default: * -- everything except SETTING files, ref. -s
  -s SETTINGS, --setting SETTINGS
                        Specify the name of SETTINGS files (yaml or json file)
                        Default: SETTINGS (case insensitive and without ext)
                        The following filename will be treat as SETTINGS file
                        SETTINGS (global settings in the folder) Start with
                        SETTINGS- (individual settings for struct) End with
                        -SETTINGS (individual settings)
  -r, --recursive       Recursive the path.
  -wf WORKFLOW, --workflow WORKFLOW
                        Specify the workflow to run. Default: get_wf_gibbs
                        (NOTE: currently, only get_wf_gibbs is supported.)
  -l, --launch          Launch the wf to launchpad
  -m [MAX_JOB], --max_job [MAX_JOB]
                        Run the job, only works when -l is specified. Default:
                        0 (Not submit to queue) 1: qlaunch singleshot (single
                        job) N(N>1): qlaunch rapidfire -m N
  -o, --write_out_wf    Write out the workflow
```

[TO TOP](#Content)