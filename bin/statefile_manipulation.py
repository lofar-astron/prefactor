#!/usr/bin/env python
import os
import argparse
import pickle
import copy

def main(statefile):
    print 'Statefile Manipulator'
    current_state = pickle.load(open(statefile, 'rb')) 
    print 'Current states:\nSetup: ', current_state[0]
    print 'List of finished steps for {0}:'.format(current_state[0]['job_name'])
    for i, item in enumerate(current_state[1]):
        name = item[0]
        try:
            name = os.path.basename(item[1]['ok.mapfile'])[:-11]
        except:
            print 'using task as name'
        print 'Step Nr.: {0}  Task: {1}  Name: {2}'.format(i+1, item[0], name)

    try:
        while(True):
            del_number_raw = raw_input('Delete last steps including number: ')
            if int(del_number_raw) < len(current_state[1])+1:
                break
            else:
                print "Not a valid step number. Check range!"
    except KeyboardInterrupt:
        print " got KeyboardInterrupt -> exiting"
        return

    # delete all steps after the given step number
    # removing a step from the middle will result in an invalid statefile
    new_state = copy.deepcopy(current_state)
    new_state[1] = new_state[1][:int(del_number_raw)-1]

    pickle.dump(new_state, open(statefile,'wb'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Delete steps from a genericpipeline-statefile")
    parser.add_argument('Statefile', help='The path to the helpfile')
    args = parser.parse_args()
    main(args.Statefile)

