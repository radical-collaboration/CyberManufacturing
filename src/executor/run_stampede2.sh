#!/bin/bash

export RADICAL_PILOT_DBURL=mongodb://chai:qwerty123@ds133816.mlab.com:33816/twoway
export RADICAL_PILOT_VERBOSE=DEBUG
export CYBER_EXECUTOR_VERBOSE=DEBUG
export RADICAL_PILOT_PROFILE=TRUE

python executor.py 2> test_stampede.log
