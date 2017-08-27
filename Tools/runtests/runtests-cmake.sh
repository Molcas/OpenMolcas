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
#               2016,2017, Ignacio Fdez. Galván                        *
#***********************************************************************
#
# runtests-cmake -- run verification tests
#
# Starting from a specified testing directory, a git clone is created for each
# file "<my-config>.cmake" which contains a cmake script that sets appropriate
# initial cache values. Then molcas is configured and built, and molcas verify
# is run with any arguments given to this script. Results are logged and mailed
# to the test page.
#
# Any arguments to this script are passed to the command 'molcas verify', so
# you can use 'runtests 044' or 'runtests caspt2'. This is to provide
# some flexibility in what is tested, whithout having to change the main
# configuration.
#
# The default testing directory is "/tmp" and the default configuration if no
# cmake files are present is to create a file "gfortran.cmake" with a default
# build type (RelWithDebInfo). This makes it easy to just run this script on any
# machine without any special arguments.
#
# Steven Vancoillie, June 2013 based on the original checkinstall scripts and
# adapted for use with git.
# Modified by Ignacio Fdez. Galván and Steven Vancoillie, including various
# improvements, parallel runs of different builds, and a cmake adaptation.
# November 2016 - March 2017: Support for several repositories.
# June 2017: Support for submodules

################################################################################
####                             CONFIGURATION                              ####
################################################################################

# directory containing:
# - files '<my-config>.cmake' containing cmake commands that set cache values,
#   the default is to create a file gfortran.flags with a line:
#     set (CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "build type")
#
# - OPTIONAL files '<my-config>.env' containing commands that should be
#   executed before running cmake (e.g. setting specific PATH or modifying
#   files). If this file is not present, this is just skipped.
#
# - OPTIONAL directories '<my-config>' which are molcas git clones (these are
#   created for each *.cmake file if the directory does not exist, so you
#   normally do not have to create these directories).
if [ -z "$TESTHOME" ]
then
    TESTHOME="/tmp"
fi

# different Molcas versions have their own branch. The main development branch
# is called "master" and every day, a branch "daily-snapshot" is available IF
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
    SUBMODULES_OPEN='External/lapack'
fi

################################################################################
####                         ADVANCED CONFIGURATION                         ####
################################################################################

# location of testpage and molcas repository
TESTPAGE='test@signe.teokem.lu.se'
UPLOADPAGE='https://molcas.altervista.org/tests/upload.php'
GITSERVER='git@git.teokem.lu.se:'
REPO='molcas-extra'
SERVER_OPEN="https://${GITLABAUTH}gitlab.com/Molcas/"
REPO_OPEN='OpenMolcas'

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
if [ -z "$CMAKE" ]
then
    CMAKE='cmake'
fi
if [ -z "$MAKE_cmd" ]
then
    MAKE_cmd='make'
    ## change the make command to build alternative targets or use different options
    #MAKE_cmd='make -j 4'
    #MAKE_cmd='make all doc'
    #MAKE_cmd='make all doc_html'
fi
if [ -z "$DRIVER" ]
then
    DRIVER='molcas'
fi
HNAME=`hostname -s`
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
# check for .cmake files with configuration flags
if ! ls *.cmake
then
    echo "no configuration files available"
    echo "creating default gfortran.cmake ..."
    echo "set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "build type")" > gfortran.cmake || exit 1
fi

checkout_clean () {
    # remove all changes to tracked files
    git reset --hard
    # quietly remove all non-tracked files
    git clean -f -x -d -q

    # over-write local branch with remote
    # first make a maintenance branch 'tmp'
    # then force-fetch the required branch
    # and check it out, then remove 'tmp'
    if git branch | grep -q "tmp"
    then
        git checkout tmp
    else
        git checkout -b tmp
    fi

    git fetch
    if git branch -r | grep -q "origin/$BRANCH"
    then
        git fetch --force origin $BRANCH:$BRANCH
        git checkout $BRANCH
        git branch -D tmp
    else
        return 1
    fi

    cd ..
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

    testconfig=${configfile%.cmake}
    if [ ! -d $testconfig ]
    then
        mkdir $testconfig
    fi

    if [ -r $testconfig.env ]
    then
        . $testconfig.env
    fi

    header=$(i=1 ; while [ $i -le $((${#testconfig}+10)) ]; do printf "%s" "#"; i=$(($i+1)); done)
    echo
    echo "$header"
    echo "#### $testconfig ####"
    echo "$header"
    echo

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
    if [ ! -d $REPO_OPEN.$BRANCH ]
    then
        ISOLD="0"
        git clone -b $BRANCH $SERVER_OPEN$REPO_OPEN $REPO_OPEN.$BRANCH || return
    fi

    for R in $REPO_OPEN $REPO
    do
        cd $R.$BRANCH || return
        if ! checkout_clean $R
        then
            echo "no branch $BRANCH available on origin ($R), skipping testing..."
            cd ../; rm -f $REPO.$BRANCH.LOCK; return
        fi
    done

    export OPENMOLCAS_DIR=`readlink -f $REPO_OPEN.$BRANCH`
    cd $REPO_OPEN.$BRANCH || return
    for sub in `git submodule foreach -q 'echo $path'` ; do
        git submodule deinit --force $sub || return 1
    done
    for sub in $SUBMODULES_OPEN ; do
        git submodule update --init --force $sub || return 1
    done

    SHA1_OPEN=`git rev-parse $BRANCH`
    MASTER_OPEN=`git rev-parse origin/master`
    echo "source tree (open) clean at $SHA1_OPEN"

    parents_open=`git log --pretty=%P -n 1 $BRANCH`
    if [ `echo $parents_open | wc -w` -eq 1 ] || [ "$SHA1_OPEN" = "$MASTER_OPEN" ]
    then
        parents_open=""
    fi

    VERSION=$(git describe --always --dirty --match "v*")
    P=$(echo $VERSION | awk -F. '{print $NF}')

    cd ../$REPO.$BRANCH || return
    for sub in `git submodule foreach -q 'echo $path'` ; do
        git submodule deinit --force $sub || return 1
    done
    for sub in $SUBMODULES ; do
        git submodule update --init --force $sub || return 1
    done

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
        if [ "$SHA1_OPEN" = "$MASTER_OPEN" ] && [ "$SHA1" = "$MASTER" ]
        then
            echo "same as master, no testing needed"
            cd ../; rm -f $REPO.$BRANCH.LOCK; return
        fi
    fi

    #### retry script ####
    ######################

    cd ../

    retry="retry.sh"
    echo "#!/bin/sh"                             >  $retry
    echo "mkdir $testconfig || exit"             >> $retry
    echo "cd $testconfig"                        >> $retry
    cat ../$testconfig.env                       >> $retry
    echo "export OPENMOLCAS_DIR=$OPENMOLCAS_DIR" >> $retry
    echo "cat << EOF > $testconfig.cmake"        >> $retry
    cat ../$testconfig.cmake                     >> $retry
    echo "EOF"                                   >> $retry
    echo "cp -r $loc/$REPO_OPEN.$BRANCH ."       >> $retry
    echo "(cd $REPO_OPEN.$BRANCH"                >> $retry
    echo "    git checkout $SHA1_OPEN"           >> $retry
    echo "    git clean -f -d -x -q)"            >> $retry
    echo "cp -r $loc/$REPO.$BRANCH ."            >> $retry
    echo "(cd $REPO.$BRANCH"                     >> $retry
    echo "    git checkout $SHA1"                >> $retry
    echo "    git clean -f -d -x -q)"            >> $retry
    echo "mkdir $REPO.$BRANCH-build"             >> $retry
    echo "cd $REPO.$BRANCH-build"                >> $retry
    echo "cmake -C ../$testconfig.cmake ../$REPO.$BRANCH > make.log 2>&1 && $MAKE_cmd VERBOSE=1 >> make.log 2>&1" >> $retry
    chmod +x $retry

    #### building ####
    ##################

    # in case of error, bark and continue with the other tests #

    echo "-- building --"

    # To ensure proper operation, delete any cache entries for an
    # existing build directory.
    if [ -d $REPO.$BRANCH-build ]
    then
        CACHEDIR=`cmake --version | awk '{print $3}'`
        rm -f $REPO.$BRANCH-build/CMakeCache.txt 2> /dev/null
        rm -f $REPO.$BRANCH-build/CMakeFiles/$CACHEDIR/* 2> /dev/null
    else
        mkdir $REPO.$BRANCH-build
    fi
    cd $REPO.$BRANCH-build || return

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
        echo "hostname    = $HNAME"     >> $outfile
        echo "uname       = $UNAME"     >> $outfile
        echo "tests       = $TESTS"     >> $outfile
        echo "Contact     = $CONTACT"   >> $outfile
        echo "Location    = $PWD"       >> $outfile
        echo "Environment:"             >> $outfile
        echo "### BEGIN ###"            >> $outfile
        cat ../../${testconfig}.env     >> $outfile
        echo "#### END ####"            >> $outfile
        echo "Initial CMake cache:"     >> $outfile
        echo "### BEGIN ###"            >> $outfile
        cat ../../${testconfig}.cmake   >> $outfile
        echo "#### END ####"            >> $outfile
    done

    # if something goes wrong it is useful to have more verbose log files
    # (change it here to allow configuring in .env files)
    MAKE_cmd="$MAKE_cmd VERBOSE=1"

    # The following flag specifies the cmake script that is executed to set
    # initial cache values.
    MY_FLAGS="-C ../../${testconfig}.cmake"

    #sed -i 's|/opt/local/bin/perl|/usr/bin/perl|' sbin/*
    # run cmake configuration, build Molcas (if you want serial output, remove
    # any "-j n" options from MAKE_cmd in the header...)
    if $CMAKE $MY_FLAGS ../$REPO.$BRANCH >> make.log 2>&1 && $MAKE_cmd >> make.log 2>&1
    then
        date >> make.log
        echo "Make - OK!" >> auto.log
    else
        echo "Make Failed! See logs." >> auto.log
        echo '************************************' >> auto.log
        echo '----------- parent details ---------' >> auto.log
        (cd ../$REPO_OPEN.$BRANCH && git checkout origin/master && update_submodules)
        for commit in $SHA1 $parents
        do
            if [ "$commit" = "$MASTER" ]
            then
                continue
            fi
            (cd ../$REPO.$BRANCH && git checkout $commit && update_submodules)
            rm -f CMakeCache.txt 2> /dev/null
            rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
            if $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1
            then
                echo ":: good $commit" >> auto.log
            else
                echo ":: fail $commit" >> auto.log
                (cd ../$REPO.$BRANCH && git log -1 --pretty=tformat:"developer: %ce" $commit) >> auto.log
            fi
        done
        (cd ../$REPO.$BRANCH && git checkout origin/master && update_submodules)
        for commit in $SHA1_OPEN $parents_open
        do
            if [ "$commit" = "$MASTER_OPEN" ]
            then
                continue
            fi
            (cd ../$REPO_OPEN.$BRANCH && git checkout $commit && update_submodules)
            rm -f CMakeCache.txt 2> /dev/null
            rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
            if $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1
            then
                echo ":: good (open) $commit" >> auto.log
            else
                echo ":: fail (open) $commit" >> auto.log
                (cd ../$REPO_OPEN.$BRANCH && git log -1 --pretty=tformat:"developer: %ce" $commit) >> auto.log
            fi
        done
        echo '************************************' >> auto.log

        # remake the original branch to save trouble for manual inspection later
        (cd ../$REPO.$BRANCH && git checkout $BRANCH && update_submodules)
        (cd ../$REPO_OPEN.$BRANCH && git checkout $BRANCH && update_submodules)
        rm -f CMakeCache.txt 2> /dev/null
        rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
        $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1

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
        (cd ../$REPO_OPEN.$BRANCH && git checkout origin/master && update_submodules)
        for commit in $SHA1 $parents
        do
            if [ "$commit" = "$MASTER" ]
            then
                continue
            fi
            (cd ../$REPO.$BRANCH && git checkout $commit && update_submodules)
            rm -f CMakeCache.txt 2> /dev/null
            rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
            if $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1 && $DRIVER verify --trap $failed_tests
            then
                echo ":: good $commit" >> auto.log
            else
                echo ":: fail $commit" >> auto.log
                (cd ../$REPO.$BRANCH && git log -1 --pretty=tformat:"developer: %ce" $commit) >> auto.log
            fi
        done
        (cd ../$REPO.$BRANCH && git checkout origin/master && update_submodules)
        for commit in $SHA1_OPEN $parents_open
        do
            if [ "$commit" = "$MASTER_OPEN" ]
            then
                continue
            fi
            (cd ../$REPO_OPEN.$BRANCH && git checkout $commit && update_submodules)
            rm -f CMakeCache.txt 2> /dev/null
            rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
            if $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1 && $DRIVER verify --trap $failed_tests
            then
                echo ":: good (open) $commit" >> auto.log
            else
                echo ":: fail (open) $commit" >> auto.log
                (cd ../$REPO_OPEN.$BRANCH && git log -1 --pretty=tformat:"developer: %ce" $commit) >> auto.log
            fi
        done
        echo '************************************' >> auto.log

        # remake the original branch to save trouble for manual inspection later
        (cd ../$REPO.$BRANCH && git checkout $BRANCH && update_submodules)
        (cd ../$REPO_OPEN.$BRANCH && git checkout $BRANCH && update_submodules)
        rm -f CMakeCache.txt 2> /dev/null
        rm -f CMakeFiles/$CACHEDIR/* 2> /dev/null
        $CMAKE $MY_FLAGS ../$REPO.$BRANCH > /dev/null 2>&1 && $MAKE_cmd > /dev/null 2>&1

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
            if ! ../$REPO_OPEN.$BRANCH/sbin/chkunprint.plx < $MESSAGE
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
for configfile in *.cmake
do
    # Run the test in a subshell, put it in the background if desired
    if [ $PARTEST -eq "1" ]
    then
        (
        logfile="$LOGDIR/$configfile.log"
        exec > "$logfile" 2>&1
        test_configfile
        )&
        echo "launched ${configfile%.cmake}, PID = $!"
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
