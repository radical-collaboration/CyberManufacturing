#!/usr/bin/env python

__copyright__ = 'Copyright 2013-2014, http://radical.rutgers.edu'
__license__   = 'MIT'
__author__    = 'Ioannis Paraskevakos'

from __future__ import division
import os
import sys
import argparse
import numpy as np
# ------------------------------------------------------------------------------
#
# READ the RADICAL-Pilot documentation: http://radicalpilot.readthedocs.org/
#
# ------------------------------------------------------------------------------

helloworld_mpi_bin  = 'helloworld_mpi.py'
helloworld_mpi_path = '%s/%s' % (os.path.abspath(os.path.dirname(__file__)),
                                 helloworld_mpi_bin)


#------------------------------------------------------------------------------
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config ',help="JSON file that contains the following information:\
        resource, cores,runtime, path to LIGGHTS,path to LIGGHTS inputs,path to PBM executable,session,\
        project.")
    parser.add_argument('--verbose', help="Verbosity level",defalut='REPORT')
    parser.add_argument('--profile', help="Profile value", defalut='FALSE')
    args = parser.parse_args()


    os.environ['RADICAL_PILOT_VERBOSE'] = args.verbose
    os.environ['RADICAL_PILOT_PROFILE'] = args.profile

    import radical.pilot as rp
    import radical.utils as ru

    # we use a reporter class for nicer output
    report = ru.LogReporter(name='radical.pilot', level=verbose)
    report.title('Getting Started (RP version %s)' % rp.version)

    # Create a new session. No need to try/except this: if session creation
    # fails, there is not much we can do anyways...
    session = rp.Session()

    # all other pilot code is now tried/excepted.  If an exception is caught, we
    # can rely on the session object to exist and be valid, and we can thus tear
    # the whole RP stack down via a 'session.close()' call in the 'finally'
    # clause...
    try:

        # read the config used for resource details
        report.info('read config')
        config = ru.read_json(args.config)
        report.ok('>>ok\n')

        report.header('submit pilots')

        # Add a Pilot Manager. Pilot managers manage one or more ComputePilots.
        pmgr = rp.PilotManager(session=session,uid=config['session'])

        # Define an [n]-core local pilot that runs for [x] minutes
        # Here we use a dict to initialize the description object
        pd_init = {
                'resource'      : config['resource'],
                'runtime'       : config['runtime'],  # pilot runtime (min)
                'exit_on_error' : True,
                'project'       : config['project'],
                'queue'         : config['queue'],
                'cores'         : config['cores']
                }
        pdesc = rp.ComputePilotDescription(pd_init)

        # Launch the pilot.
        pilot = pmgr.submit_pilots(pdesc)


        report.header('submit units')

        # Register the ComputePilot in a UnitManager object.
        umgr = rp.UnitManager(session=session)
        umgr.add_pilots(pilot)

        # Create a workload of ComputeUnits. 
        # Each compute unit runs a MPI test application.
        report.info('create DEM unit description\n\t')


        # create a new CU description, and fill it.
        # Here we don't use dict initialization.
        cud = rp.ComputeUnitDescription()
        cud.environment    = ['PATH='+config['pathtoLIGGHTS']+':$PATH']
        cud.pre_exec       = ['mkdir CSVs','mkdir post','mkdir restart']
        cud.executable     = 'lmp_micstam'
        cud.arguments      = ['-in','in.*' ]
        cud.input_staging  = []
        cud.output_staging = [{'source': 'unit:///post/collision200000.atom',
                               'target': 'pilot:///collision200000.atom',
                               'mode'  : rp.LINK},
                              {'source': 'unit:///post/impact200000.atom',
                               'target': 'pilot:///impact200000.atom',
                               'mode'  : rp.LINK}]
        cud.cores          = config['cores']
        cud.mpi            = True
        cuds.append(cud)
        report.progress()
        report.ok('>>ok\n')

        # Submit the previously created ComputeUnit descriptions to the
        # PilotManager. This will trigger the selected scheduler to start
        # assigning ComputeUnits to the ComputePilots.
        units = umgr.submit_units(cuds)

        # Wait for all compute units to reach a final state (DONE, CANCELED or FAILED).
        report.header('gather results')
        umgr.wait_units()
    
        report.info('\n')
        collision = {'source': 'pilot:///post/collision200000.atom',
                     'target': 'unit:///sampledumpfiles/collision200000.%d_200'%config['cores'],
                     'mode'  : rp.LINK}

        impact = {'source': 'pilot:///impact200000.atom',
                  'target': 'unit:///sampledumpfiles/impact200000.%d_200'%config['cores'],
                  'mode'  : rp.LINK}
        # create a new CU description, and fill it.
        # Here we don't use dict initialization.
        cud = rp.ComputeUnitDescription()
        cud.environment    = ['PATH='+config['pathtoPBM']+':$PATH']
        cud.pre_exec       = ['mkdir sampledumpfiles']
        cud.executable     = 'model.out'
        cud.input_staging  = [collision,impact]
        cud.cores          = 16
        cud.mpi            = True
        cuds.append(cud)
        report.progress()
        report.ok('>>ok\n')

    except Exception as e:
        # Something unexpected happened in the pilot code above
        report.error('caught Exception: %s\n' % e)
        raise

    except (KeyboardInterrupt, SystemExit) as e:
        # the callback called sys.exit(), and we can here catch the
        # corresponding KeyboardInterrupt exception for shutdown.  We also catch
        # SystemExit (which gets raised if the main threads exits for some other
        # reason).
        report.warn('exit requested\n')

    finally:
        # always clean up the session, no matter if we caught an exception or
        # not.  This will kill all remaining pilots.
        report.header('finalize')
        session.close()

    report.header()


#-------------------------------------------------------------------------------

