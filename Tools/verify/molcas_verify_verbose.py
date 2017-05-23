#!/usr/bin/env python
#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) Thomas Bondo Pedersen                                  *
#***********************************************************************
#
# Author: Thomas Bondo Pedersen
#         Centre for Theoretical and Computational Chemistry
#         University of Oslo, Norway
#

import os, sys
from optparse import OptionParser, OptionGroup

def getLines(filename):
    '''Return lines of file filename as a list, each item containing a line
       without the newline character'''
    f=open(filename,'r')
    lines=[line.split('\n')[0] for line in f.readlines()]
    f.close()
    return lines

def findOccurrences(item,lines,case_sensitive=False):
    '''Return indices of lines containing item'''
    if case_sensitive:
        return [i for i in range(len(lines)) if item in lines[i]]
    else:
        return [i for i in range(len(lines)) if item.lower() in lines[i].lower()]

def isHappyLanding(lines):
    last=lines[findOccurrences('/rc=',lines)[-1]]
    rc=last.split('/rc=')[1].split()[0]
    try:
        rc=int(rc)
    except ValueError:
        rc=1
    return rc==0

# set up option parser
usage='''%prog [options] log-files'''
description='''Parse molcas verify log-files and report failed and, optionally, skipped test results. Typical usage, assuming %prog is located in your PATH,
with molcas verification: molcas verify -- %prog test*.log'''
parser=OptionParser(usage=usage,description=description)
group=OptionGroup(parser,'Output control options')
group.add_option('--no-color',
                 action='store_false',
                 dest='use_color',
                 default=True,
                 help='Turn off color printing (red=failed, yellow=passed with skipped failures, green=passed)')
group.add_option('--skipped','--report-skipped',
                 action='store_true',
                 dest='report_skipped',
                 default=False,
                 help='Report skipped tests')
parser.add_option_group(group)

# parse options
(options,args)=parser.parse_args()

# set colors on or off
dsh='----------'
if not options.use_color:
    RED=''
    YELLOW=''
    GREEN=''
    ENDC=''
else:
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    GREEN='\033[0;32m'
    ENDC='\033[0m'

# Look for these words in log file to find failed/skipped tests
Failed='Failed'
Skipped='Skipped'
# process one log file at a time
for fil in args:
    lines=getLines(fil)  # read log file and get list of lines without newline characters
    if len(lines)<1:
        print RED+os.path.split(fil)[1].split('.')[0]+ENDC
        sys.stdout.flush()
    else:
        # find line numbers where modules end (excluding check and auto)
        calculation_failed=not isHappyLanding(lines)
        ids=[i for i in findOccurrences('--- Stop Module:',lines,case_sensitive=True) if 'check' not in lines[i] and 'auto' not in lines[i]]
        # get sections of the log file where the check is done for each module
        data={}
        sections=[]
        range_ids=range(len(ids))
        ids+=[-1]
        for k in range_ids:
            i=ids[k]
            raw_section=lines[i].split()[3]
            count=-1
            Stop=False
            while not Stop:
                count+=1
                section=raw_section+' '+str(count)
                if not data.has_key(section):
                    Stop=True
                    j=ids[k+1]
                    data[section]=lines[i:j]
                    sections.append(section)
        # get failed and skipped tests
        n_failed=0
        n_skipped=0
        report_failed={}
        report_skipped={}
        for section in sections:
            report_failed[section]=[RED+data[section][i]+ENDC for i in findOccurrences(Failed,data[section],case_sensitive=True)]
            report_skipped[section]=[YELLOW+data[section][i]+ENDC for i in findOccurrences(Skipped,data[section],case_sensitive=True)]
            n_failed+=len(report_failed[section])
            n_skipped+=len(report_skipped[section])
        # print result
        if n_failed>0:
            print RED+dsh+os.path.split(fil)[1].split('.')[0]+dsh+ENDC
            for section in sections:
                if len(report_failed[section])>0:
                    print RED+' '+section.split()[0]+ENDC
                    for txt in report_failed[section]:
                        print '  '+txt
                    if options.report_skipped:
                        for txt in report_skipped[section]:
                            print '  '+txt
                else:
                    if len(report_skipped[section])>0:
                        print YELLOW+' '+section.split()[0]+ENDC
                        if options.report_skipped:
                            for txt in report_skipped[section]:
                                print '  '+txt
                    else:
                        print GREEN+' '+section.split()[0]+ENDC
        else:
            if n_skipped>0:
                if calculation_failed:
                    print RED+dsh+os.path.split(fil)[1].split('.')[0]+dsh+ENDC
                else:
                    print YELLOW+dsh+os.path.split(fil)[1].split('.')[0]+dsh+ENDC
                for section in sections:
                    if len(report_skipped[section])>0:
                        print YELLOW+' '+section.split()[0]+ENDC
                        if options.report_skipped:
                            for txt in report_skipped[section]:
                                print '  '+txt
                    else:
                        print GREEN+' '+section.split()[0]+ENDC
            else:
                if calculation_failed:
                    print RED+dsh+os.path.split(fil)[1].split('.')[0]+dsh+ENDC
                else:
                    print GREEN+dsh+os.path.split(fil)[1].split('.')[0]+dsh+ENDC
                for section in sections:
                    print GREEN+' '+section.split()[0]+ENDC
        sys.stdout.flush()
