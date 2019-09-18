#!/bin/sh
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
# Copyright (C) 2013, Steven Vancoillie                                *
#               2016-2019, Ignacio Fdez. Galván                        *
#***********************************************************************
#
# runtests-configure -- run verification tests
#
# Starting from a specified testing directory, a git clone is created for each
# file "<my-config>.flag" which contains lines with configure flags. Then
# molcas is built with those flags and molcas verify is run with any arguments
# give to this script. Results are logged and mailed to the test page.
#
# Any arguments to this script are passed to the command 'molcas verify', so
# you can use 'runtests 044' or 'runtests caspt2'. This is to provide
# some flexibility in what is tested, whithout having to change the main
# configuration.
#
# The default testing directory is "/tmp" and the default configuration if no
# flag files are present is to create a file "gfortran.flags" with no special
# configuration flags (just "-serial"). This makes it easy to just run this
# script on any machine without any special arguments.
#
# Steven Vancoillie, June 2013 based on the original checkinstall scripts ands
# adapted for use with git.
# Modified by Ignacio Fdez. Galván, December 2016 - March 2017: Support for
# several repositories.
# June 2017: Support for submodules
# August 2017: Use fetch
# July-August 2018: Remove OpenMolcas testing, this should now only handle
#   molcas-extra

################################################################################
####                             CONFIGURATION                              ####
################################################################################

# directory containing:
# - files '<my-config>.flags' containing configure flags if there are no files,
#   the default is to create a file gfortran.flags with a line "-serial
#   -packages" (configuration always includes flag -noprompt)
#
# - OPTIONAL files '<my-config>.cmd' containing commands that should be
#   executed before running configure (e.g. setting specific PATH or modifying
#   files). If this file is not present, this is just skipped.
#
# - OPTIONAL files '<my-config>.rte' that should replace molcas.rte before
#   running any tests. If this file is not present, this is just skipped.
#
# - OPTIONAL directories '<my-config>' which are molcas git clones (these are
#   created for each *.flags file if the directory does not exist, so you
#   normally do not have to create these directories).
if [ -z "$TESTHOME" ]
then
    TESTHOME="/tmp"
fi

# different Molcas versions have their own branch. The main development branch
# is called "master" and every day, a branch "master-testing" is available IF
# there were new personal branches marked for testing. The default job of a
# test setup is to run daily tests on this branch. You can also test "master"
# itself, or if there are stable versions, they will have branches "rel-8.0" and
# "rel-8.0-testing" respectively (names to be determined).
if [ -z "$BRANCH" ]
then
    BRANCH="daily-snapshot"
fi

# mail program, if mail not available you can try perlmail
if [ -z "$MAIL_cmd" ]
then
    #MAIL_cmd='perlmail'
    MAIL_cmd="mail"
fi

# contact information (YOUR name and email address)
if [ -z "$CONTACT" ]
then
    CONTACT='Firstname Lastname youremail@domain'
fi

# you can set a global PATH here (e.g. if run through cron)
#PATH=''

# if you want to get notification by mail - add it into RECIPIENT
if [ -z "$RECIPIENT" ]
then
    RECIPIENT=''
fi

# submodules to be updated by default (space-sparated list)
if [ -z "$SUBMODULES" ]
then
    SUBMODULES=''
fi
if [ -z "$SUBMODULES_OPEN" ]
then
    SUBMODULES_OPEN=''
fi

################################################################################
####                         ADVANCED CONFIGURATION                         ####
################################################################################

# location of testpage and molcas repository
TESTPAGE='test@signe.teokem.lu.se'
UPLOADPAGE='https://molcas.altervista.org/tests/upload.php'
GITSERVER='git@git.teokem.lu.se:'
REPO='molcas-extra'
REPO_OPEN='openmolcas'

RECIPIENT="$RECIPIENT $TESTPAGE"

## set this to 1 if you want all the different configurations to be run simultaneously
PARTEST=0

## some commands (could be different on your platform)
#GUNZIP='gzip -d'
#UNTAR='tar -xf'
if [ -z "$FOLD" ]
then
    FOLD='fold -w 70'
fi
if [ -z "$MAKE_cmd" ]
then
    MAKE_cmd='make'
fi
if [ -z "$DRIVER" ]
then
    DRIVER='molcas'
fi
HNAME=`hostname -s`
FHNAME=`hostname -A`
UNAME=`uname -a`
DATE=`date +%F_%T`
LANG=C

################################################################################

# you can specify tests to run on the command
# line or change it here
TESTS="$*"
if [ -z "$TESTS" ]
then
    TESTS=".all"
fi

# sanity check on the test directory
if [ ! -d $TESTHOME ]
then
    mkdir $TESTHOME || exit 1
fi

# redirect output to logfile
LOGDIR="$TESTHOME/checklog/"
if [ ! -d $LOGDIR ]
then
    mkdir $LOGDIR || exit 1
fi
logfile="$LOGDIR/$DATE"
exec > "$logfile" 2>&1

# before exiting, mail the logfile to the user running the tests
#trap "$MAIL_cmd -s \"Molcas testing results from $HNAME\" $CONTACT" EXIT

cd $TESTHOME || exit 1

# existing configuration flag files?
# check for .flags files with configuration flags
if ! ls *.flags
then
    echo "no configuration files available"
    echo "creating default gfortran.flags ..."
    echo "-serial" > gfortran.flags || exit 1
fi

checkout_clean () {
    # remove all changes to tracked files
    git reset --hard || return -1
    # quietly remove all non-tracked files
    git clean -f -x -d -q || return -1

    # over-write local branch with remote
    # first make a maintenance branch 'tmp'
    # then force-fetch the required branch
    # and check it out, then remove 'tmp'
    if git branch | grep -q "tmp"
    then
        git checkout tmp || return -1
    else
        git checkout -b tmp || return -1
    fi

    git fetch || return -1
    if git branch -r | grep -q "origin/$1"
    then
        git fetch --force origin $1:$1 || return -1
        git checkout $1 || return -1
        git branch -D tmp || return -1
    else
        return 1
    fi
}

update_submodules () {
    for sub in `git submodule foreach -q 'echo $path'` ; do
        git submodule update --init --force $sub || return 1
    done
}

cut_log () {
    if [ `wc -l < "$1"` -gt 30000 ] ; then
        head -n 15000 "$1" > "$1.cut"
        echo "~~~ file too long, lines removed ~~~" >> "$1.cut"
        tail -n 15000 "$1" >> "$1.cut"
        mv "$1" "$1.long"
        mv "$1.cut" "$1"
    fi
}

test_configfile () {
    if [ ! -r "$configfile" ]
    then
        echo
        echo "error: cannot read $configfile, skipping test..."
        return
    fi

    testconfig=${configfile%.flags}
    if [ ! -d $testconfig ]
    then
        mkdir $testconfig
    fi

    if [ -r $testconfig.cmd ]
    then
        . $testconfig.cmd
    fi

    header=$(i=1 ; while [ $i -le $((${#testconfig}+10)) ]; do printf "%s" "#"; i=$(($i+1)); done)
    echo
    echo "$header"
    echo "#### $testconfig ####"
    echo "$header"
    echo

    # get configuration flags
    MY_FLAGS="-noprompt"
    for line in `cat $configfile`
    do
        MY_FLAGS="$MY_FLAGS $line"
    done

    cd $testconfig

    #### preparation ####
    #####################

    echo "<= fetching branch $BRANCH =>"; echo

    if [ -f $REPO.$BRANCH.LOCK ]
    then
        echo "$BRANCH is still busy, skipping ..."
        return
    else
        touch $REPO.$BRANCH.LOCK
        trap "rm -f $PWD/$REPO.$BRANCH.LOCK" EXIT
    fi

    ISOLD="1"
    if [ ! -d $REPO.$BRANCH ]
    then
        ISOLD="0"
        git clone -b $BRANCH $GITSERVER$REPO $REPO.$BRANCH || return
    fi
    if [ -d $REPO_OPEN ]
    then
        rm -rf $REPO_OPEN || return
    fi

    cd $REPO.$BRANCH || return
    checkout_clean $BRANCH
    rc=$?
    if [ $rc -lt 0 ]
    then
        echo "error checking out $BRANCH from $REPO, skipping testing..."
        cd ../; rm -f $REPO.$BRANCH.LOCK; return
    elif [ $rc -gt 0 ]
    then
        echo "no branch $BRANCH available on origin ($REPO), skipping testing..."
        cd ../; rm -f $REPO.$BRANCH.LOCK; return
    fi

    for sub in `git submodule foreach -q 'echo $path'` ; do
        git submodule deinit --force $sub || return 1
    done
    for sub in $SUBMODULES ; do
        git submodule update --init --force $sub || return 1
    done

    ./fetch OPENMOLCAS
    cd ../$REPO_OPEN || return
    for sub in `git submodule foreach -q 'echo $path'` ; do
        git submodule deinit --force $sub || return 1
    done
    for sub in $SUBMODULES_OPEN ; do
        git submodule update --init --force $sub || return 1
    done
    export OPENMOLCAS_DIR=`readlink -f $PWD`
    SHA1_OPEN=`git rev-parse HEAD`
    VERSION=$(git describe --always --dirty --match "v*")
    P=$(echo $VERSION | awk -F. '{print $NF}')
    cd ../$REPO.$BRANCH || return
    rm -f fetch.log

    SHA1=`git rev-parse $BRANCH`
    MASTER=`git rev-parse origin/master`
    echo "source tree clean at $SHA1"

    # get parents of this commit for testing in case of failure
    # since git rev-list is only in recent versions, use git log
    # if there is a single parent, or the branch coincides with master,
    # ignore any parent
    parents=`git log --pretty=%P -n 1 $BRANCH`
    if [ `echo $parents | wc -w` -eq 1 ] || [ "$SHA1" = "$MASTER" ]
    then
        parents=""
    fi

    # get a human readable version for use with testpage
    # we do this here first in case configure/make fail
    VERSION=$(git describe --always --dirty --match "v*")
    V=$(echo $VERSION | awk -F. '{print $1"."$2}' | tr -d v)
    P="$P & "$(echo $VERSION | awk -F. '{print $NF}')

    # if the directories had been cloned (and presumably tested) before,
    # test the branches only if they differ from the master (no mail sent)
    if [ $ISOLD -eq "1" ]
    then
        if [ "$SHA1" = "$MASTER" ]
        then
            echo "same as master, no testing needed"
            cd ../; rm -f $REPO.$BRANCH.LOCK; return
        fi
    fi

    #### retry script ####
    ######################

    cd ../

    retry="retry.sh"
    echo "#!/bin/sh"                  >  $retry
    echo "mkdir $testconfig || exit"  >> $retry
    echo "cd $testconfig"             >> $retry
    cat ../$testconfig.cmd            >> $retry
    echo "cp -r $PWD/$REPO.$BRANCH ." >> $retry
    echo "cd $REPO.$BRANCH"           >> $retry
    echo "git checkout $SHA1"         >> $retry
    echo "git clean -f -d -x -q"      >> $retry
    echo "./configure $MY_FLAGS > make.log 2>&1 && $MAKE_cmd >> make.log 2>&1" >> $retry
    chmod +x $retry

    #### building ####
    ##################

    # in case of error, bark and continue with the other tests #

    echo "-- building --"

    cd $REPO.$BRANCH

    # for automatic reporting, taken from checkinstall
    #   make.log: contains output of build
    #   auto.log: contains output of tests
    printf "=:=: makelog: %s _ %s\n" \
        $HNAME $testconfig > make.log
    printf "=:=: autolog: %s _ %s\n" \
        $HNAME $testconfig > auto.log

    # some header information that goes in both files:
    for outfile in make.log auto.log
    do
        echo "MOLCAS version $V patch level $P" >> $outfile
        date >> $outfile
        echo "SHA1        = $SHA1"      >> $outfile
        echo "SHA1 (open) = $SHA1_OPEN" >> $outfile
        echo "hostname    = $FHNAME"    >> $outfile
        echo "uname       = $UNAME"     >> $outfile
        echo "tests       = $TESTS"     >> $outfile
        echo "Contact     = $CONTACT"   >> $outfile
        echo "Location    = $PWD"       >> $outfile
        echo "Flags       = $MY_FLAGS"  >> $outfile
        echo "Environment:"             >> $outfile
        echo "### BEGIN ###"            >> $outfile
        cat ../../${testconfig}.cmd     >> $outfile
        echo "#### END ####"            >> $outfile
        echo "Flags file:"              >> $outfile
        echo "### BEGIN ###"            >> $outfile
        cat ../../${testconfig}.flags   >> $outfile
        echo "#### END ####"            >> $outfile
        echo "runtests version: 190918" >> $outfile
    done

    # run configure and make
    if ./configure $MY_FLAGS >> make.log 2>&1 && $MAKE_cmd >> make.log 2>&1
    then
        date >> make.log
        echo "Make - OK!" >> auto.log
    else
        echo "Make Failed! See logs." >> auto.log
        echo '************************************' >> auto.log
        echo '----------- parent details ---------' >> auto.log
        for commit in $SHA1 $parents
        do
            if [ "$commit" = "$MASTER" ]
            then
                continue
            fi
            (git reset --hard && git checkout $commit && update_submodules && ./fetch OPENMOLCAS && cd ../$REPO_OPEN && update_submodules)
            make distclean >/dev/null 2>&1
            if ./configure $MY_FLAGS >/dev/null 2>&1 && $MAKE_cmd >/dev/null 2>&1
            then
                echo ":: good $commit" >> auto.log
            else
                echo ":: fail $commit" >> auto.log
                git log -1 --pretty=tformat:"developer: %ce" $commit >> auto.log
            fi
        done
        echo '************************************' >> auto.log

        # remake the original branch to save trouble for manual inspection later
        (git reset --hard && git checkout $BRANCH && update_submodules && ./fetch OPENMOLCAS && cd ../$REPO_OPEN && update_submodules)
        make distclean >/dev/null 2>&1
        ./configure $MY_FLAGS >/dev/null 2>&1 && $MAKE_cmd >/dev/null 2>&1

        # end logfiles with the date
        date >> make.log
        date >> auto.log

        # final line to indicate end to Mr. Robot
        echo "=:=:" >> make.log
        echo "=:=:" >> auto.log

        # send logfiles to Molcas testpage
        for RCPT in $RECIPIENT
        do
            for MESSAGE in make.log auto.log
            do
                cut_log $MESSAGE
                curl --form "file=@$MESSAGE" $UPLOADPAGE || curl --ciphers ecdhe_ecdsa_aes_256_sha --form "file=@$MESSAGE" $UPLOADPAGE
                echo "sending mail"
                $FOLD -s $MESSAGE | $MAIL_cmd -s $DATE $RCPT
                mv $MESSAGE.long $MESSAGE 2> /dev/null
            done
        done

        # print last part of the build and skip remainder of test
        tail make.log
        cd ../; rm -f $REPO.$BRANCH.LOCK; return
    fi

    echo '----------- Input ---------' >> auto.log

    if [ -r ../../$testconfig.rte ]
    then
        cp ../../$testconfig.rte molcas.rte
    fi

    #### testing ####
    #################

    # in case of error, bark and continue with the other tests #

    echo "-- testing --"

    if $DRIVER verify $TESTS
    then
        cat test/result >> auto.log
        # there might still be failures for 8xx series, log them
        for failed in `ls test/failed/*.out 2> /dev/null`
        do
            echo '************************************' >> auto.log
            echo '--- Test' ${failed##*/} >> auto.log
            tail -100 "$failed" >> auto.log
        done
    else
        cat test/result >> auto.log
        # log the last part of the tests
        for failed in `ls test/failed/*.out 2> /dev/null`
        do
            echo '************************************' >> auto.log
            echo '--- Test' ${failed##*/} >> auto.log
            tail -100 "$failed" >> auto.log
        done

        grep 'Failed!' test/result
        failed_tests=`$DRIVER verify --failed --flatlist`

        echo '************************************' >> auto.log
        echo '----------- parent details ---------' >> auto.log
        for commit in $SHA1 $parents
        do
            if [ "$commit" = "$MASTER" ]
            then
                continue
            fi
            (git reset --hard && git checkout $commit && update_submodules && ./fetch OPENMOLCAS && cd ../$REPO_OPEN && update_submodules)
            make distclean >/dev/null 2>&1
            if ./configure $MY_FLAGS >/dev/null 2>&1 && $MAKE_cmd >/dev/null 2>&1 && $DRIVER verify --trap $failed_tests
            then
                echo ":: good $commit" >> auto.log
            else
                echo ":: fail $commit" >> auto.log
                git log -1 --pretty=tformat:"developer: %ce" $commit >> auto.log
            fi
        done
        echo '************************************' >> auto.log

        # remake the original branch to save trouble for manual inspection later
        (git reset --hard && git checkout $BRANCH && update_submodules && ./fetch OPENMOLCAS && cd ../$REPO_OPEN && update_submodules)
        make distclean >/dev/null 2>&1
        ./configure $MY_FLAGS >/dev/null 2>&1 && $MAKE_cmd >/dev/null 2>&1
    fi

    # end logfiles with the date
    date >> make.log
    date >> auto.log

    # final line to indicate end to Mr. Robot
    echo "=:=:" >> make.log
    echo "=:=:" >> auto.log

    # send logfiles to Molcas testpage
    for RCPT in $RECIPIENT
    do
        for MESSAGE in make.log auto.log
        do
            cut_log $MESSAGE
            curl --form "file=@$MESSAGE" $UPLOADPAGE || curl --ciphers ecdhe_ecdsa_aes_256_sha --form "file=@$MESSAGE" $UPLOADPAGE
            # remove garbage from any output
            if ! ../$REPO_OPEN/sbin/chkunprint.plx < $MESSAGE
            then
                cp $MESSAGE ${MESSAGE}.gbg
                tr -c '\11\12\15\40-\176' ';' < ${MESSAGE}.gbg > $MESSAGE
            fi
            # send mail to testpage
            echo "sending mail"
            $FOLD -s $MESSAGE | $MAIL_cmd -s $DATE $RCPT
            mv $MESSAGE.long $MESSAGE 2> /dev/null
        done
    done

    cd ../
    rm -f $REPO.$BRANCH.LOCK
}

# go through all tests in sequence
files=""
for configfile in *.flags
do
    # Run the test in a subshell, put it in the background if desired
    if [ $PARTEST -eq "1" ]
    then
        (
        logfile="$LOGDIR/$configfile.log"
        exec > "$logfile" 2>&1
        test_configfile
        )&
        echo "launched ${configfile%.flags}, PID = $!"
    else
        ( test_configfile )
    fi

    files="$files $configfile"

# end loop over config files
done

# If the tests were run in the background, wait for them to
# finish and collect their output
if [ $PARTEST -eq "1" ]
then
    echo "Done! Subshells are now running the verification..."

    wait

    for configfile in $files
    do
        logfile="$LOGDIR/$configfile.log"
        cat "$logfile"
        rm -f "$logfile"
    done
fi

exit 0
