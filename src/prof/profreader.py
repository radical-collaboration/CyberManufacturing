import pandas as pd
import numpy as np
import argparse

def profile_reader(profile=None):

    prof = open(profile).readlines()
    profDF = pd.DataFrame(columns=['uid','start','end'])
    for profLine in prof:
        if ('unit' in profLine) and \
           ('update_request' in profLine) and \
           ('AGENT_EXECUTING' in profLine) and \
           ('AGENT_EXECUTING_PENDING' not in profLine):
           uid = profLine.split(',')[2]
           unit_start = float(profLine.split(',')[0])
           profDF.loc[len(profDF)]=[uid,unit_start,None]
        elif ('unit' in profLine) and \
             ('update_request' in profLine) and \
             ('AGENT_STAGING_OUTPUT_PENDING' in profLine):
           uid = profLine.split(',')[2]
           unit_end = float(profLine.split(',')[0])
           indx = profDF[profDF['uid']==uid].index
           profDF.loc[indx[0]] = ['unit.000000',profDF.get_value(index=indx[0],col='start'),unit_end]

    return profDF

if __name__ == '__main__':    

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename',help='Name of the profile file')
    args = parser.parse_args()

    df = profile_reader(profile=args.filename)

    start = df['start'].values[0]
    end   = df['end'].values[-1]
    total_time = end-start
    overhead = end-start
    for idx, values in df.iterrows():
        print 'unit:',idx,'duration:',values['end']-values['start']
        overhead -= values['end']-values['start']

    print 'Total Time:',total_time,'Overhead:',overhead