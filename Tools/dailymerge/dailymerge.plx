#!/usr/bin/perl
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
# dailymerge.plx:
#
# script that is run daily to:
# 1) check if daily-snapshot differs from master, if so
#    then check if verification has succeeded:
#    - yes: advance master to daily-snapshot, tag as stable
#    - no: reset daily-snapshot to master, tag failed branches
# 2) merge all main nickname branches into daily-snapshot
#    except those that are labeled 'failed-' or contain changes
#    to restricted directories (e.g. sbin/)
#
# Steven Vancoillie, Lund, May-June 2013
#
# started as perl incarnation of dailymerge.sh for better pattern matching, and
# changed a lot since (checking test page, handling sbin/, automatic selection
# of branches in case of merge conflict, manually exclude branches).
#
# Ignacio Fdez. Galván, Uppsala, November 2016 - April 2017
#
# modified to handle two repositories, and other improvements

# Perl modules
use warnings;
use strict;
use Getopt::Std;
use List::Util qw(first);

# set hot pipes, flushes output
$| = 1;

# option parsing
my %opt = ();
getopts("hdliufsg", \%opt);

my $help = $opt{h};
my $is_local = $opt{l};
my $is_interactive = $opt{i};
my $update_only = $opt{u};
my $forced_update = $opt{f};
my $snapshot_only = $opt{s};
my $goodriddens = $opt{g};

sub usage {
    print <<USAGE;

 Usage: dailymerge.plx [OPTION]...

        Program to check the results of automated verification,
        update the master branch if there are no failures, and
        make a new snapshot if there are new commits.

    -i     interactive: ask permission before executing any git
           commands that could alter the repository

    -l     local: do not push anything back to origin

    -u     update only:

           check if tests passed, advance the master branch

    -f     force an update (used with -u)

           treat 'DOWN' as 'FAIL', resulting in an update
           regardless of the state of the testing machines.

    -s     snapshot only:

           make a new daily snapshot

    -g     update as usual, but in case of failure pass on the
           separate verified (good) branches to the snapshot

    -h     help: print usage information

USAGE
    exit 0;
}

# subroutine to handle git commands
sub git {
    if ( $is_local && ($_[0] eq 'push') ) {
        print "(LOCAL, NO PUSH)> git @_\n";
        return 0;
    } else {
        print "> git @_\n";
        my $doit = "y";
        if ( $is_interactive ) {
          print "ok? (y/n)\n";
          chomp ($doit = <STDIN>);
        }
        system "git", @_ if ( lc $doit eq 'y' );
    }
}

################################################################################
####                             CONFIGURATION                              ####
################################################################################

# test-page access
my $username='**********@**********';
my $password='**********';
my $testpage='http://molcas.org/dev/test';

# location of the local repository in which the commits will be applied
my $localrepo = q(/home/gitupdater/molcas-extra);
my $localrepo_open = q(/home/gitupdater/openmolcas);

# location of the directory containing the source code snapshots
my $snapshots = q(/home/gitupdater/snapshots);

# directory receiving the log files:
my $logroot = q(/home/gitupdater/log);

# email address for sending log files:
my $usermail = q(**********@**********);

# branches authorized to change sbin/ go here
my @admin_devs = qw(
                       valera
                       stevenv
                       ignacio
                  );

################################################################################
####                             GENERAL SETUP                              ####
################################################################################

&usage if $help;

# set date
chomp(my $date = `date +%y%m%d-%H%M`);
my $current_time = time;
my @now = localtime;

### openmolcas ###

## change current directory to the local repository
chdir $localrepo_open or die 'Error: cannot change to repo (open) directory';

print "DAILY MERGE\n";
print "===========\n";
print "Fetch latest updates from origin (open):\n";
&git("fetch", "-p");

# get current hashes for master and daily-snapshot
chomp(my $master_open = `git rev-parse origin/master`);
chomp(my $daily_open  = `git rev-parse origin/daily-snapshot`);
print " origin/master          = $master_open\n";
print " origin/daily-snapshot  = $daily_open\n";

# Check if the daily snapshot is an actual merge (multiple parents) or if
# it corresponds to someone's branch. This could make a difference for sending
# out emails.
my @parents_open = split (/ +/, `git log --pretty=%P -n 1 $daily_open`);
my $number_of_parents_open = @parents_open;

### molcas-extra ###

## change current directory to the local repository
chdir $localrepo or die 'Error: cannot change to repo (extra)directory';

print "Fetch latest updates from origin (extra):\n";
&git("fetch", "-p");

# get current hashes for master and daily-snapshot
chomp(my $master = `git rev-parse origin/master`);
chomp(my $daily  = `git rev-parse origin/daily-snapshot`);
print " origin/master          = $master\n";
print " origin/daily-snapshot  = $daily\n";

# Check if the daily snapshot is an actual merge (multiple parents) or if
# it corresponds to someone's branch. This could make a difference for sending
# out emails.
my @parents = split (/ +/, `git log --pretty=%P -n 1 $daily`);
my $number_of_parents = @parents;

goto SNAPSHOT if $snapshot_only;

################################################################################
####                        UPDATE MASTER BRANCHES                          ####
################################################################################

UPDATE:
print "\n";
print "Update master branches\n";
print "----------------------\n";

foreach ($localrepo, $localrepo_open) {
    chdir $_ or die "Error: cannot change to directory $_";
    &git("checkout", "master", "-q");
    &git("merge", "--ff-only", "origin/master", "-q");
}

# global variable to hold the merge heads, in order to allow these to be set
# during the update section in case we have to make a snapshot with verified
# branches (-g option)
my @merge_heads;
my @merge_heads_open;

$number_of_parents = 0 if ($master eq $daily);
$number_of_parents_open = 0 if ($master_open eq $daily_open);

my %failed_configs;
my %failed_parents;
my %failed_parents_open;
unless (($master eq $daily) and ($master_open eq $daily_open)) {

    # Get hashes of failed snapshots from certified servers.

    # First, get the list of certified servers
    my @certified = ();
    my $cert_file = '/home/gitupdater/certificate';
    open CERT, '<', $cert_file or die "Error: failed to open certificate file\n";
    while (<CERT>) {
        next if /^#/;
        chomp;
        push @certified, $_;
    }
    close CERT;

    # For each server in the list, fetch the auto.log file and first look for
    # the SHA1 and make sure it matches the daily snapshot. If it does, all is
    # well and the server has completed verification. Then search at the end of
    # the file for the lines that hold the hashes of the parents that failed
    # verification.

    print "Checking certified servers:\n";
    my $sleepy_count = 0;
    my $digest = "/tmp/digest-$date";
    open DIGEST, '>', $digest or die "Error: failed to open file\n";

    foreach my $config (@certified) {
        # start section in digest file
        print DIGEST "<==== Tester: $config ====>\n";

        # get the auto.log file from the server
        printf " %-40s %3s ", $config, "...";
        my $status = system 'wget', '-N', '-q', '-P', '/tmp',
                   "--user=$username",
                   "--password=$password",
                   "$testpage/$config.txt";
        die "Error: failed to get file: $!\n" unless $status == 0;

        my $autolog = "/tmp/$config.txt";
        open TESTPAGE, '<', $autolog or die "Error: failed to open file\n";

        my $sha1;
        my $sha1_open;
        my $down = 0;
        my $contact;
        while (<TESTPAGE>) {
            next unless /^SHA1  /;
            chomp;
            $sha1 = (split)[2];
            while (<TESTPAGE>) {
                next unless /^SHA1 \(open\)/;
                chomp;
                $sha1_open = (split)[3];
                last;
            }
            if (($sha1 ne $daily) or ($sha1_open ne $daily_open)) {
                while (<TESTPAGE>) {
                    next unless /^Contact/;
                    chomp;
                    ($contact) = ((split(/=/))[1] =~ m/.*?<?(\S*@\S*)>?/);
                    last;
                }
                print "DOWN\n";
                $down = 1;
                $sleepy_count++;
                # treat this as a failure if we don't want to
                # postpone the update (with the -f flag)
                $failed_configs{$config}++ if $forced_update;
                my $msg = "/tmp/msg";
                open MSG, '>', $msg or die "Error: failed to open file\n";
                print MSG "The machine $config has not reported the test results\n";
                print MSG "for the latest daily snapshots (extra: $daily, open: $daily_open).\n";
                print MSG "Please check the machine and fix the problem as soon\n";
                print MSG "as possible.\n\n";
                print MSG "Common causes are:\n";
                print MSG "- Verification was slow and didn't finish on time, or stalled\n";
                print MSG "- A previous unfinished verification left a .LOCK file\n";
                print MSG "- A problem occurred sending the results by mail\n";
                close MSG;
                system("mail -s \"[molcas-git] certified machine did not finish\" $contact < $msg");
                print "Mail sent to $contact\n";
            }
            last;
        }

        die "Error: no SHA1 hash in auto.log\n" unless $sha1;
        die "Error: no SHA1 (open) hash in auto.log\n" unless $sha1_open;

        unless ($down) {
            my $badtests = "/tmp/badtests";
            open BADTESTS, '>', $badtests or die "Error: failed to open file\n";
            print BADTESTS "Today's snapshot as well as your separate branch did not pass\n";
            print BADTESTS "verification on the $config platform.\n\n";
            print BADTESTS "Below is a list of _all_ failures on this platform, not just yours.\n";
            print BADTESTS "- for make failures, check the right link on the test page for more info.\n";
            print BADTESTS "- for both make and verification failures, the left link on the test page\n";
            print BADTESTS "  provides a tail of either the make log or failed test logs.\n\n";
            print BADTESTS "Summary of failures:\n\n";

            # print lines of failed tests
            my $failed_verification = 0;
            while (<TESTPAGE>) {
                last if /^\*+$/;
                print DIGEST if /Failed!/;
                print BADTESTS if /Failed!/;
                $failed_verification++ if /^make.*failed/i;
                $failed_verification++ if /^standard.*failed/i;
                $failed_verification++ if /^extra.*failed/i;
            }

            # print tails of the job outputs
            while (<TESTPAGE>) {
                last if /parent details/;
                print DIGEST;
            }

            my @failed = ();
            my @failed_open = ();
            while (<TESTPAGE>) {
                next unless s/^:: fail\s+//;
                chomp;
                if (s/\(open\)\s+//) {
                    push @failed_open, $_;
                } else {
                    push @failed, $_;
                }
            }

            close BADTESTS;

            # printing
            if ($failed_verification or @failed or @failed_open) {
                $failed_configs{$config}++;
                print "FAIL\n";
                if ((($number_of_parents > 1) and not @failed) and
                    (($number_of_parents_open > 1) and not @failed_open)) {
                    print STDERR "WARNING!! snapshot failed, but parents are OK!!\n" unless (@failed or @failed_open);
                }
            } else {
                print "PASS\n";
            }

            # emailing
            # The distinction between actual merges (multiple parents) and single parents
            # does not make much sense when every "snapshot" is a merge of two repositories
            # so just send a mail to every committer whose branch failed
            chdir $localrepo or die 'Error: cannot change to repo (extra) directory';
            if ($number_of_parents > 0) {
                foreach (@failed) {
                    $failed_parents{$_}++;
                    chomp(my $committer = `git show -s --pretty=format:%ce $_`);
                    print "  bad parent (extra) $_ ($committer)\n";
                    system("mail -s \"[molcas] your snapshot (extra) failed on $config\" $committer < $badtests") unless ($committer eq 'gitupdater@signe.teokem.lu.se');
                }
            }
            chdir $localrepo_open or die 'Error: cannot change to repo (open) directory';
            if ($number_of_parents_open > 0) {
                foreach (@failed_open) {
                    $failed_parents_open{$_}++;
                    chomp(my $committer = `git show -s --pretty=format:%ce $_`);
                    print "  bad parent (open) $_ ($committer)\n";
                    system("mail -s \"[molcas] your snapshot (open) failed on $config\" $committer < $badtests") unless ($committer eq 'gitupdater@signe.teokem.lu.se');
                }
            }

            unlink $badtests;
        }

        close TESTPAGE;

        unlink $autolog;
    }
    close DIGEST;

    if ($sleepy_count > 0) {
        print STDERR "Some servers have not completed testing.\n";
        unless ($forced_update or %failed_configs) {
            print STDERR "Please check this, now exiting dailymerge...\n";
            exit 1;
        }
    }
} else {
    print "master (extra & open) already at daily snapshot, no update needed\n";
    goto SNAPSHOT;
}

unless ($master eq $daily) {

    chdir $localrepo or die 'Error: cannot change to repo (extra) directory';

    # check if the tag we will use is still available
    my $tag_master = "master-x$date";
    die "Error: tag $tag_master already exists\n" if `git tag -l $tag_master`;

    if (%failed_configs) {
        print STDERR "Failed configurations:\n";
        print STDERR " $_\n" foreach keys %failed_configs;
    } else {
        print "Verification complete\n";
    }

    # Now, unless there are any failed snapshots, fast-forward the
    # master branch to the daily snapshot. If not, tag the failed
    # branches as such and reset the daily branch back to master.

    unless (%failed_configs) {
        # forward master to daily snapshot
        &git ("merge", "--ff-only", "origin/daily-snapshot");
        &git ("tag", "-a", $tag_master, "-m", "verified snapshot");

        # push changes to remote origin
        print qq(Pushing tags and updating master on remote (extra):\n);
        &git("push", "--tags");
        &git("push", "origin", "master:master");

        # update link to stable snapshot
        chomp (my $version = `git tag -l --points-at master v*`);
        system ("(cd $snapshots && ln -sf molcas-extra-$version.tar.gz molcas-extra-dev-stable.tgz)");

    } else {
        # go through SHA1 of failed commits and tag them as such
        my $i = 0;
        # the merge commit is listed as a "parent"
        #&git("tag", "FAILED-x$date.".$i++, $daily);
        foreach (keys %failed_parents) {
            if ($_ eq $daily) {
                die "Error: tag FAILED-x$date already exists\n" if `git tag -l FAILED-x$date`;
                &git("tag", "FAILED-x$date", $_);
            } else {
                $i++;
                die "Error: tag FAILED-x$date.$i already exists\n" if `git tag -l FAILED-x$date.$i`;
                &git("tag", "FAILED-x$date.$i", $_);
            }
        }

        # reset daily-snapshot
        &git("checkout", "daily-snapshot");
        &git("reset", "--hard", "master");
        &git("checkout", "master");

        # push changes to remote origin
        print qq(Pushing tags and resetting daily-snapshot on remote (extra):\n);
        &git("push", "--tags");
        &git("push", "-f", "origin", "daily-snapshot:daily-snapshot");

        # If we need to make a snapshot from verified branches (-g option), then
        # set the list of merge heads here by getting the parents of the failed
        # snapshot and exclude those that have failed.
        if ($goodriddens) {
            if ($number_of_parents > 1) {
                my @parents = split (/ +/, `git log --pretty=%P -n 1 $daily`);
                foreach (@parents) {
                    chomp; push @merge_heads, $_ unless $failed_parents{$_};
                }
            # If there is only one parent, check if it actually failed
            } else {
                push @merge_heads, $daily unless $failed_parents{$daily};
            }
        }
    }
} else {
    print "master (extra) already at daily snapshot, no update needed\n";
}

unless ($master_open eq $daily_open) {

    chdir $localrepo_open or die 'Error: cannot change to repo (open) directory';

    # check if the tag we will use is still available
    my $tag_master = "master-o$date";
    die "Error: tag $tag_master already exists\n" if `git tag -l $tag_master`;

    if (%failed_configs) {
        print STDERR "Failed configurations:\n";
        print STDERR " $_\n" foreach keys %failed_configs;
    } else {
        print "Verification complete\n";
    }

    # Now, unless there are any failed snapshots, fast-forward the
    # master branch to the daily snapshot. If not, tag the failed
    # branches as such and reset the daily branch back to master.

    unless (%failed_configs) {
        # forward master to daily snapshot
        &git ("merge", "--ff-only", "origin/daily-snapshot");
        &git ("tag", "-a", $tag_master, "-m", "verified snapshot");

        # push changes to remote origin
        print qq(Pushing tags and updating master on remote (open):\n);
        &git("push", "--tags");
        &git("push", "origin", "master:master");

        # update link to stable snapshot
        chomp (my $version = `git tag -l --points-at master v*`);
        system ("(cd $snapshots && ln -sf openmolcas-$version.tar.gz openmolcas-dev-stable.tgz)");

    } else {
        # go through SHA1 of failed commits and tag them as such
        my $i = 0;
        # the merge commit is listed as a "parent"
        #&git("tag", "FAILED-o$date.".$i++, $daily_open);
        foreach (keys %failed_parents_open) {
            if ($_ eq $daily_open) {
                die "Error: tag FAILED-o$date already exists\n" if `git tag -l FAILED-o$date`;
                &git("tag", "FAILED-o$date", $_);
            } else {
                $i++;
                die "Error: tag FAILED-o$date.$i already exists\n" if `git tag -l FAILED-o$date.$i`;
                &git("tag", "FAILED-o$date.$i", $_);
            }
        }

        # reset daily-snapshot
        &git("checkout", "daily-snapshot");
        &git("reset", "--hard", "master");
        &git("checkout", "master");

        # push changes to remote origin
        print qq(Pushing tags and resetting daily-snapshot on remote (open):\n);
        &git("push", "--tags");
        &git("push", "-f", "origin", "daily-snapshot:daily-snapshot");

        # If we need to make a snapshot from verified branches (-g option), then
        # set the list of merge heads here by getting the parents of the failed
        # snapshot and exclude those that have failed.
        if ($goodriddens) {
            if ($number_of_parents_open > 1) {
                my @parents_open = split (/ +/, `git log --pretty=%P -n 1 $daily_open`);
                foreach (@parents_open) {
                    chomp; push @merge_heads_open, $_ unless $failed_parents_open{$_};
                }
            # If there is only one parent, check if it actually failed
            } else {
                push @merge_heads_open, $daily_open unless $failed_parents_open{$daily_open};
            }
        }
    }
} else {
    print "master (open) already at daily snapshot, no update needed\n";
}

# now checkout daily-snapshot
SNAPSHOT:

exit 0 if $update_only;

################################################################################
####                         CREATE DAILY SNAPSHOT                          ####
################################################################################

# At this point, local master and daily-snapshot HAVE to have the same hash,
# otherwise something is wrong.  Either master was moved to daily-snapshot
# because verification was OK, or daily-snapshot was reset to master because
# verification failed. Either way, they must have the same SHA1 hash.

### molcas-extra ###

chdir $localrepo or die 'Error: cannot change to repo (extra) directory';

print "\n";
print "Checkout daily-snapshot branch (extra)\n";
print "--------------------------------------\n";

&git("checkout", "daily-snapshot");

chomp($master = `git rev-parse master`);
chomp($daily  = `git rev-parse daily-snapshot`);
print "master          = $master\n";
print "daily-snapshot  = $daily\n";

if ($master ne $daily) {
    print 'Error: master and daily-snapshot do not match';

    exit 1;
}

### openmolcas ###

chdir $localrepo_open or die 'Error: cannot change to repo (open) directory';

print "\n";
print "Checkout daily-snapshot branch (open)\n";
print "-------------------------------------\n";

&git("checkout", "daily-snapshot");

chomp($master_open = `git rev-parse master`);
chomp($daily_open  = `git rev-parse daily-snapshot`);
print "master          = $master_open\n";
print "daily-snapshot  = $daily_open\n";

if ($master_open ne $daily_open) {
    print 'Error: master and daily-snapshot (open) do not match';

    exit 1;
}

### molcas-extra ###

chdir $localrepo or die 'Error: cannot change to repo (extra) directory';

print "\n";
print "Update daily-snapshot branch (extra)\n";
print "------------------------------------\n";

my @unauthorized_branches;
my @failed_branches;

unless ($goodriddens) {
    # Get a list of all remote branches that have not been merged yet and
    # remove HEAD, master, daily-snapshot or any branches having a tag
    # beginning with 'FAILED-'. For all valid branches, remove any branches
    # that do not belong to an admin and introduce changes in sbin/, because
    # this is not allowed.

    my @remote_branches = `git branch --remote --no-merged`;
    foreach (@remote_branches) {
        chomp; s/^\s+//; s/\s+$//;
    }

    my @valid_branches;
    my %exclude_tip;
    foreach my $branch (@remote_branches) {
        (my $branchname = $branch) =~ s!^origin/!!;

        # record branches to be excluded
        (my $nickname = $branchname) =~ s/-.*$//;
        if ($branchname eq $nickname."-EXCLUDE") {
            chomp(my $sha1 = `git rev-parse $branch`);
            $exclude_tip{$nickname} = $sha1;
        } else {
            $exclude_tip{$nickname} = "";
        }

        # skip nickname-... branches
        next if $branchname =~ m/\w+-/;

        # skip malformed names (like "origin/HEAD -> origin/master")
        next if $branchname =~ m/\s/;

        # remove failed branches
        my @tags = `git tag --points-at $branch`;
        if (grep /^FAILED-/, @tags) {
            push @failed_branches, $branch;
            next;
        }

        # if there is the special "patchmaster" branch, then take only that
        # branch, assign it to the merge heads and proceed immediately to the
        # merge phase.
        if ($branchname =~ m/patchmaster/) {
            @merge_heads = ($branch);
            goto PATCHMASTER;
        }

        # manually skip/include branches
        #next if $branchname =~ m/devname/;
        #next unless $branchname =~ m/devname/;

        push @valid_branches, $branch;
    }

    foreach my $branch (@valid_branches) {
        (my $nickname = $branch) =~ s!^origin/!!;

        # skip nickname-EXCLUDE branches
        chomp(my $sha1 = `git rev-parse $branch`);
        if ($exclude_tip{$nickname} eq $sha1) {
            print "EXCLUDE $branch\n";
            next;
        }

        # check for unauthorized changes
        chomp(my $merge_base = `git merge-base HEAD $branch`);
        my @files = `git diff --name-only $merge_base $branch`;
        unless (grep /^$nickname$/, @admin_devs) {
            if ( grep /^sbin\//, @files ) {
                print "Warning: $nickname not allowed to change files\n";
                print "         in sbin/, excluding branch $branch...\n";
                push @unauthorized_branches, $branch;
            } else {
                push @merge_heads, $branch;
            }
        } else {
            push @merge_heads, $branch;
        }
    }
} else {
    # The merge heads were inherited from the update section, but be paranoid
    # and check again that none of them have a FAILED- tag pointing at them!
    foreach my $commit (@merge_heads) {
        my @tags = `git tag --points-at $commit`;
        if (grep /^FAILED-/, @tags) {
            print "Internal Error: inherited merge head was supposed to be OK\n";
            print "                offending SHA1 (extra) hash = $commit\n";
            exit 1;
        }
    }
}

# manually skip emails (for testing purposes)
#goto PATCHMASTER;

# Send mail for unauthorized (changes in sbin) branches
my $msg = "/tmp/msg";
open MSG, '>', $msg or die "Error: failed to open file\n";
print MSG "Today, your branch was not included in our daily snapshot (extra),\n";
print MSG "because it contained changes in the sbin/ directory.\n";
print MSG "This is restricted and only allowed for administrators.\n\n";
print MSG "If you wish to apply the changes, please contact one of the\n";
print MSG "the administrators.\n\n";
print MSG "You can undo your changes in sbin/ by running:\n";
print MSG "  git checkout origin/master -- sbin\n";
print MSG "and committing the result as usual.\n";
close MSG;
foreach my $branch (@unauthorized_branches) {
    chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
    system("mail -s \"[molcas-git] there were unauthorized changes\" $committer < $msg");
    #system("mail -s \"[molcas-git] there were unauthorized changes for $committer\" $usermail < $msg");
    print "Mail sent to $committer\n";
}

# Check and report old unmerged branches, only on Mondays PM
if (($now[6] == 1) && ($now[2] >= 12)) {
    open MSG, '>', $msg or die "Error: failed to open file\n";
    print MSG "Your branch in molcas-extra failed two weeks ago or earlier. It has not\n";
    print MSG "been updated since then, so it is not being considered for merges with\n";
    print MSG "the master.  If you want your changes included, please merge with the\n";
    print MSG "current master branch, fix any problem that caused your branch to fail,\n";
    print MSG "and push your branch to the repository. If you don't do anything, you\n";
    print MSG "will continue to get this message every week, you can avoid this by\n";
    print MSG "creating a branch called <nickname>-EXCLUDE, where <nickname> is the\n";
    print MSG "name of your normal branch.\n";
    close MSG;
    print "\nExcluded (FAILED) branches:\n";
    foreach my $branch (@failed_branches) {
        chomp(my $sha1 = `git rev-parse $branch`);
        chomp(my $ts = `git show -s --format=%ct $branch`);
        my $age = int(($current_time-$ts)/86400);
        print "$branch: age: $age days old\n";
        if ($age >= 14) {
            chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
            system("mail -s \"[molcas-git] you have old unmerged changes\" $committer < $msg");
            #system("mail -s \"[molcas-git] $committer has old unmerged changes\" $usermail < $msg");
            print "Mail sent to $committer\n";
        }
    }
}

PATCHMASTER:
my $number_of_merge_attempts = 0;

MERGE:
$number_of_merge_attempts++;

if (@merge_heads) {
    print "branches/commits selected for merging (extra):\n";
    print " $_\n" foreach @merge_heads;
} else {
    print "No valid unmerged branches (extra), done.\n";
    goto MERGE_DONE;
}

# Merge the selected branches into daily-snapshot at master.
#
# Now that we have the needed branches in @merge_heads, merge them into the
# daily-snapshot branch, at the tip of master. We can just merge everything as
# git itself will make sure nothing unnecessary is merged (already merged
# branches are ignored, so we do not need to check for that!).

my $status = &git("merge", "-q", "-m", "Merge developer branches into daily-snapshot", "--log", @merge_heads);
if ($status == 0) {
    # merging went fine

    chomp(my $merged = `git rev-parse HEAD`);
    goto MERGE_DONE if $merged eq $daily;

    # if a version tag already exist here, we don't create a new one
    my @tags = `git tag --points-at $merged`;
    chomp(my $tag_match = first { $_ =~ m/^v8\.1\.x/ } @tags) if @tags;
    if ($tag_match) {

        # link old snapshot
        system ("(cd $snapshots && ln -sf molcas-extra-$tag_match.tar.gz molcas-extra-dev-daily.tgz)");

    } else {

        # first, check if certain tags exist that shouldn't
        my $tag_daily  = "v8.3.x$date";
        die "Error: tag $tag_daily already exists\n" if `git tag -l $tag_daily`;

        &git("tag", "-a", $tag_daily, "-m", "daily snapshot for testing");

        # prepare snapshot
        system ("(cd $localrepo_open ; MOLCAS_SOURCE=$localrepo molcas export && mv $localrepo/molcas-$tag_daily.tar.gz $snapshots/molcas-extra-$tag_daily.tar.gz)");
        system ("(cd $snapshots && ln -sf molcas-extra-$tag_daily.tar.gz molcas-extra-dev-daily.tgz)");

    }

    # apply the local changes now to the remote
    print "Pushing tags and daily-snapshot to remote origin (extra):\n";
    &git("push", "--tags");
    &git("push", "origin", "daily-snapshot:daily-snapshot");

} else {
    # there was a problem with merging

    my $signal = $status & 127;
    my $retval = $status >> 8;
    print STDERR "Error: merge failed (extra) with return code $retval\n";

    # get conflicting files
    my @conflicts = `git diff --name-only --cached --diff-filter=U`;

    # identify merge heads that touch the conflicting files
    my %bad_branches;
    foreach my $branch (@merge_heads) {
        chomp(my $merge_base = `git merge-base HEAD $branch`);
        my @files = `git diff --name-only $merge_base $branch`;
        foreach my $conflict (@conflicts) {
            chomp $conflict;
            if (grep /\Q$conflict\E/, @files) {
                push @{ $bad_branches{$branch} }, $conflict;
            }
        }
    }

    # if the number of bad branches is 1, then assume it is in conflict with
    # the master branch or it is the only remaining branch of the bunch. It
    # should just be deleted later, so let it stay in the bad_branches hash.
    if (keys %bad_branches == 1) {
        print "single bad branch: ", keys %bad_branches, "\n";
    } else {
        # for multiple bad branches, take the first one that has merged the latest
        # origin/master branch. If none are found, traverse all master snapshots till
        # there is one that has been merged by a branch, and then take that branch
        my @master_tags = `git tag --list master-* | sort -r`;
      SELECT_GOOD_BRANCH:
        foreach my $master_tag (@master_tags) {
            chomp $master_tag;
            chomp (my $master_hash = `git rev-parse $master_tag^{commit}`);
            foreach my $branch (keys %bad_branches) {
                chomp(my $merge_base = `git merge-base $branch $master_hash`);
                if ($master_hash eq $merge_base) {
                    # remove one of the chosen branches from the bad branches
                    delete $bad_branches{$branch};
                    last SELECT_GOOD_BRANCH;
                }
                ;
            }
        }
    }

    # manually skip emails (for testing purposes)
    #goto RESET;

    my $msg = "/tmp/msg";
    open MSG, '>', $msg or die "Error: failed to open file\n";
    print MSG "Today, your branch was not included in our daily snapshot (extra),\n";
    print MSG "because it conflicted with origin/master or another branch.\n\n";
    print MSG "This is the list of conflicting branches/files:\n";
    foreach my $branch (keys %bad_branches) {
        chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
        print MSG "\nconflicting files in: $branch ($committer)\n";
        foreach my $conflict (@{ $bad_branches{$branch} }) {
            print MSG " $conflict\n";
        }
    }
    print MSG "\n";
    print MSG "In order to resolve this issue, you should merge the\n";
    print MSG "origin/master branch tomorrow and resolve any conflicts.\n";
    print MSG "You can find how to do this on our developer's DokuWiki.\n\n";
    print MSG "If you have no time to do this, don't worry, you will\n";
    print MSG "just keep receiving this message but otherwise no-one\n";
    print MSG "else is affected by this issue.\n\n";
    print MSG "Please note that conflicts can be avoided by checking\n";
    print MSG "if a file has been changed by someone but was not yet\n";
    print MSG "merged. From your current branch, execute:\n\n";
    print MSG "git log --no-merges --remotes ^HEAD -- path/to/file\n";
    close MSG;
    foreach my $branch (keys %bad_branches) {
        chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
        system("mail -s \"[molcas-git] there were conflicts during merge\" $committer < $msg");
        #system("mail -s \"[molcas-git] there were conflicts during merge for $committer\" $usermail < $msg");
        print "Mail sent to $committer\n";
    }

  RESET:
    &git("reset", "--hard");

    # finally, remove the bad branches from the merge heads and then
    # jump back to the merge procedure, unless we have already tried
    # it a couple of times (set to 2 tries, increase if necessary).
    my @old_heads = @merge_heads;
    @merge_heads = ();
    foreach my $head (@old_heads) {
        push @merge_heads, $head unless $bad_branches{$head};
    }
    if ($number_of_merge_attempts > 2) {
        print "-> I reached the maximum number of merge attempts (extra).\n";
        print "-> Either increase the number of attempts or wait\n";
        print "-> until some developers have manually resolved the\n";
        print "-> merge conflicts in their branch.\n";
    } else {
        print "-> trying new merge\n";
        goto MERGE;
    }

}

MERGE_DONE:

### openmolcas ###

chdir $localrepo_open or die 'Error: cannot change to repo (open) directory';

print "\n";
print "Update daily-snapshot branch (open)\n";
print "-----------------------------------\n";

my @unauthorized_branches_open;
my @failed_branches_open;

unless ($goodriddens) {
    # Get a list of all remote branches that have not been merged yet and
    # remove HEAD, master, daily-snapshot or any branches having a tag
    # beginning with 'FAILED-'. For all valid branches, remove any branches
    # that do not belong to an admin and introduce changes in sbin/, because
    # this is not allowed.

    my @remote_branches = `git branch --remote --no-merged`;
    foreach (@remote_branches) {
        chomp; s/^\s+//; s/\s+$//;
    }

    my @valid_branches;
    my %exclude_tip;
    foreach my $branch (@remote_branches) {
        (my $branchname = $branch) =~ s!^origin/!!;

        # record branches to be excluded
        (my $nickname = $branchname) =~ s/-.*$//;
        if ($branchname eq $nickname."-EXCLUDE") {
            chomp(my $sha1 = `git rev-parse $branch`);
            $exclude_tip{$nickname} = $sha1;
        } else {
            $exclude_tip{$nickname} = "";
        }

        # skip nickname-... branches
        next if $branchname =~ m/\w+-/;

        # skip malformed names (like "origin/HEAD -> origin/master")
        next if $branchname =~ m/\s/;

        # remove failed branches
        my @tags = `git tag --points-at $branch`;
        if (grep /^FAILED-/, @tags) {
            push @failed_branches_open, $branch;
            next;
        }

        # if there is the special "patchmaster" branch, then take only that
        # branch, assign it to the merge heads and proceed immediately to the
        # merge phase.
        if ($branchname =~ m/patchmaster/) {
            @merge_heads_open = ($branch);
            goto PATCHMASTER_open;
        }

        # manually skip/include branches
        #next if $branchname =~ m/devname/;
        #next unless $branchname =~ m/devname/;

        push @valid_branches, $branch;
    }

    foreach my $branch (@valid_branches) {
        (my $nickname = $branch) =~ s!^origin/!!;

        # skip nickname-EXCLUDE branches
        chomp(my $sha1 = `git rev-parse $branch`);
        if ($exclude_tip{$nickname} eq $sha1) {
            print "EXCLUDE $branch\n";
            next;
        }

        # check for unauthorized changes
        chomp(my $merge_base = `git merge-base HEAD $branch`);
        my @files = `git diff --name-only $merge_base $branch`;
        unless (grep /^$nickname$/, @admin_devs) {
            if ( grep /^sbin\//, @files ) {
                print "Warning: $nickname not allowed to change files\n";
                print "         in sbin/, excluding branch $branch...\n";
                push @unauthorized_branches_open, $branch;
            } else {
                push @merge_heads_open, $branch;
            }
        } else {
            push @merge_heads_open, $branch;
        }
    }
} else {
    # The merge heads were inherited from the update section, but be paranoid
    # and check again that none of them have a FAILED- tag pointing at them!
    foreach my $commit (@merge_heads_open) {
        my @tags = `git tag --points-at $commit`;
        if (grep /^FAILED-/, @tags) {
            print "Internal Error: inherited merge head was supposed to be OK\n";
            print "                offending SHA1 (open) hash = $commit\n";
            exit 1;
        }
    }
}

# manually skip emails (for testing purposes)
#goto PATCHMASTER_open;

# Send mail for unauthorized (changes in sbin) branches
$msg = "/tmp/msg";
open MSG, '>', $msg or die "Error: failed to open file\n";
print MSG "Today, your branch was not included in our daily snapshot (open),\n";
print MSG "because it contained changes in the sbin/ directory.\n";
print MSG "This is restricted and only allowed for administrators.\n\n";
print MSG "If you wish to apply the changes, please contact one of the\n";
print MSG "the administrators.\n\n";
print MSG "You can undo your changes in sbin/ by running:\n";
print MSG "  git checkout origin/master -- sbin\n";
print MSG "and committing the result as usual.\n";
close MSG;
foreach my $branch (@unauthorized_branches_open) {
    chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
    system("mail -s \"[molcas-git] there were unauthorized changes\" $committer < $msg");
    #system("mail -s \"[molcas-git] there were unauthorized changes for $committer\" $usermail < $msg");
    print "Mail sent to $committer\n";
}

# Check and report old unmerged branches, only on Mondays PM
if (($now[6] == 1) && ($now[2] >= 12)) {
    open MSG, '>', $msg or die "Error: failed to open file\n";
    print MSG "Your branch in openmolcas failed two weeks ago or earlier. It has not\n";
    print MSG "been updated since then, so it is not being considered for merges with\n";
    print MSG "the master.  If you want your changes included, please merge with the\n";
    print MSG "current master branch, fix any problem that caused your branch to fail,\n";
    print MSG "and push your branch to the repository. If you don't do anything, you\n";
    print MSG "will continue to get this message every week, you can avoid this by\n";
    print MSG "creating a branch called <nickname>-EXCLUDE, where <nickname> is the\n";
    print MSG "name of your normal branch.\n";
    close MSG;
    print "\nExcluded (FAILED) branches:\n";
    foreach my $branch (@failed_branches_open) {
        chomp(my $ts = `git show -s --format=%ct $branch`);
        my $age = int(($current_time-$ts)/86400);
        print "$branch: age: $age days old\n";
        if ($age >= 14) {
            chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
            system("mail -s \"[molcas-git] you have old unmerged changes\" $committer < $msg");
            #system("mail -s \"[molcas-git] $committer has old unmerged changes\" $usermail < $msg");
            print "Mail sent to $committer\n";
        }
    }
}

PATCHMASTER_open:
$number_of_merge_attempts = 0;

MERGE_open:
$number_of_merge_attempts++;

if (@merge_heads_open) {
    print "branches/commits selected for merging (open):\n";
    print " $_\n" foreach @merge_heads_open;
} else {
    print "No valid unmerged branches (open), done.\n";
    goto MERGE_DONE_open;
}

# Merge the selected branches into daily-snapshot at master.
#
# Now that we have the needed branches in @merge_heads_open, merge them into the
# daily-snapshot branch, at the tip of master. We can just merge everything as
# git itself will make sure nothing unnecessary is merged (already merged
# branches are ignored, so we do not need to check for that!).

$status = &git("merge", "-q", "-m", "Merge developer branches into daily-snapshot", "--log", @merge_heads_open);
if ($status == 0) {
    # merging went fine

    chomp(my $merged = `git rev-parse HEAD`);
    goto MERGE_DONE_open if $merged eq $daily_open;

    # if a version tag already exist here, we don't create a new one
    my @tags = `git tag --points-at $merged`;
    chomp(my $tag_match = first { $_ =~ m/^v8\.1\.o/ } @tags) if @tags;
    if ($tag_match) {

        # link old snapshot
        system ("(cd $snapshots && ln -sf openmolcas-$tag_match.tar.gz openmolcas-dev-daily.tgz)");

    } else {

        # first, check if certain tags exist that shouldn't
        my $tag_daily  = "v8.3.o$date";
        die "Error: tag $tag_daily already exists\n" if `git tag -l $tag_daily`;

        &git("tag", "-a", $tag_daily, "-m", "daily snapshot for testing");

        # prepare snapshot
        system ("(cd $localrepo_open ; MOLCAS_SOURCE=$localrepo_open molcas export && mv $localrepo_open/molcas-$tag_daily.tar.gz $snapshots/openmolcas-$tag_daily.tar.gz)");
        system ("(cd $snapshots && ln -sf openmolcas-$tag_daily.tar.gz openmolcas-dev-daily.tgz)");

    }

    # apply the local changes now to the remote
    print "Pushing tags and daily-snapshot to remote origin (open):\n";
    &git("push", "--tags");
    &git("push", "origin", "daily-snapshot:daily-snapshot");

} else {
    # there was a problem with merging

    my $signal = $status & 127;
    my $retval = $status >> 8;
    print STDERR "Error: merge failed with return code $retval\n";

    # get conflicting files
    my @conflicts = `git diff --name-only --cached --diff-filter=U`;

    # identify merge heads that touch the conflicting files
    my %bad_branches;
    foreach my $branch (@merge_heads_open) {
        chomp(my $merge_base = `git merge-base HEAD $branch`);
        my @files = `git diff --name-only $merge_base $branch`;
        foreach my $conflict (@conflicts) {
            chomp $conflict;
            if (grep /\Q$conflict\E/, @files) {
                push @{ $bad_branches{$branch} }, $conflict;
            }
        }
    }

    # if the number of bad branches is 1, then assume it is in conflict with
    # the master branch or it is the only remaining branch of the bunch. It
    # should just be deleted later, so let it stay in the bad_branches hash.
    if (keys %bad_branches == 1) {
        print "single bad branch: ", keys %bad_branches, "\n";
    } else {
        # for multiple bad branches, take the first one that has merged the latest
        # origin/master branch. If none are found, traverse all master snapshots till
        # there is one that has been merged by a branch, and then take that branch
        my @master_tags = `git tag --list master-* | sort -r`;
      SELECT_GOOD_BRANCH:
        foreach my $master_tag (@master_tags) {
            chomp $master_tag;
            chomp (my $master_hash = `git rev-parse $master_tag^{commit}`);
            foreach my $branch (keys %bad_branches) {
                chomp(my $merge_base = `git merge-base $branch $master_hash`);
                if ($master_hash eq $merge_base) {
                    # remove one of the chosen branches from the bad branches
                    delete $bad_branches{$branch};
                    last SELECT_GOOD_BRANCH;
                }
                ;
            }
        }
    }

    # manually skip emails (for testing purposes)
    #goto RESET_open;

    my $msg = "/tmp/msg";
    open MSG, '>', $msg or die "Error: failed to open file\n";
    print MSG "Today, your branch was not included in our daily snapshot (open),\n";
    print MSG "because it conflicted with origin/master or another branch.\n\n";
    print MSG "This is the list of conflicting branches/files:\n";
    foreach my $branch (keys %bad_branches) {
        chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
        print MSG "\nconflicting files in: $branch ($committer)\n";
        foreach my $conflict (@{ $bad_branches{$branch} }) {
            print MSG " $conflict\n";
        }
    }
    print MSG "\n";
    print MSG "In order to resolve this issue, you should merge the\n";
    print MSG "origin/master branch tomorrow and resolve any conflicts.\n";
    print MSG "You can find how to do this on our developer's DokuWiki.\n\n";
    print MSG "If you have no time to do this, don't worry, you will\n";
    print MSG "just keep receiving this message but otherwise no-one\n";
    print MSG "else is affected by this issue.\n\n";
    print MSG "Please note that conflicts can be avoided by checking\n";
    print MSG "if a file has been changed by someone but was not yet\n";
    print MSG "merged. From your current branch, execute:\n\n";
    print MSG "git log --no-merges --remotes ^HEAD -- path/to/file\n";
    close MSG;
    foreach my $branch (keys %bad_branches) {
        chomp(my $committer = `git show -s --pretty=format:%ce $branch`);
        system("mail -s \"[molcas-git] there were conflicts during merge\" $committer < $msg");
        #system("mail -s \"[molcas-git] there were conflicts during merge for $committer\" $usermail < $msg");
        print "Mail sent to $committer\n";
    }

  RESET_open:
    &git("reset", "--hard");

    # finally, remove the bad branches from the merge heads and then
    # jump back to the merge procedure, unless we have already tried
    # it a couple of times (set to 2 tries, increase if necessary).
    my @old_heads = @merge_heads_open;
    @merge_heads_open = ();
    foreach my $head (@old_heads) {
        push @merge_heads_open, $head unless $bad_branches{$head};
    }
    if ($number_of_merge_attempts > 2) {
        print "-> I reached the maximum number of merge attempts (open).\n";
        print "-> Either increase the number of attempts or wait\n";
        print "-> until some developers have manually resolved the\n";
        print "-> merge conflicts in their branch.\n";
    } else {
        print "-> trying new merge\n";
        goto MERGE_open;
    }

}

MERGE_DONE_open:

exit 0;
