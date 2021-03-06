# Introduction

This README file contains basic instructions on how to execute the RADICAL-Pilot 
script that couples the DEM and PBM execution.

It is suggested to read [RADICAL-Pilot's documentation](https://radicalpilot.readthedocs.io/en/latest/).
The following instructions are more of quick setup and run.

## RADICAL-Pilot Installation
If you are GCC Python execute the following:

```
vitrualenv rp
source rp/bin/activate
pip install radical.pilot
pip install numpy
```

In case of Anaconda Python:

```
conda create -y -p rp python=2.7 numpy
source activate rp
pip install radical.pilot
```

## Execution
Before executing the script, please create a json file that will hold a number 
of information to configure the execution. Those are:

*    resource: RADICAL-Pilot Resource label,
*    cores: Number of cores that will be used by the agent,
*    diameter: The diameter of the DEM particles
*    runtime: Total Requested runtime,
*    queue: The queue that will be used,
*    pathtoLIGGGHTS: Path to LIGGGHTS executable on the remote resource,
*    pathtoLIGGHTSinputs: Path to DEM imput folder,
*    pathtoPBMexecutable: Path to PBM executable,
*    session: A meaningfull name (it can also be null, RADICAL-Pilot will autoname it),
*    project: Allocation number

An example would look like the json file in this folder.

Create a MongoDB database on MLab or some other MongoService and do:

```
export RADICAL_PILOT_DBURL=<mongodburl>
```


When finished run as
```
python rp_script.py config.json
```

For additional information about the scripts arguments execute:
```
python rp_script.py -h
````

When the script has finished executing, the user can inspect their output data in RADICAL-Pilot's sandbox. The sandbox 
exist in the selected remote resource, under the folder ```radical.pilot.sandbox```.  A folder with the session name will exist.
For example, if the session was named ```foo``` in the json config file, the result of the execution will be under ```radical.pilot.sandbox/foo```. DEM results are under ```unit.000000``` and PBM results under ```unit.000001``` folder.
