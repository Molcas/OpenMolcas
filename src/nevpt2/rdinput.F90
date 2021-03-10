!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! MOLCAS wrapper for Celestino Angeli's NEVPT2 code
! rdinput reads input from the MOLCAS input file
!
!*******************

subroutine rdinput(refwfnfile)

! use global variables directly from the NEVPT2 program
use nevpt2_cfg, only: igelo, MultGroup, no_pc, nr_frozen_orb, nr_states, rdm_distributed, rdm_path, rdm_read, skip_effective_ham, &
                      skip_koopro_molcas
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(out) :: refwfnfile
character(len=180) :: Line, key
character(len=9001) :: dline
character(len=9001) :: frozen_str
integer(kind=iwp) :: LuSpool, iError, i, isplit
integer(kind=iwp), external :: isFreeUnit
logical(kind=iwp), external :: next_non_comment
character(len=180), external :: Get_Ln

! Initial values

refwfnfile = 'JOBIPH'

LuSpool = isFreeUnit(18)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'NEVPT2')
iError = -1

!> set it to one in order to make sure NEVPT2 runs also with just the basic input
!> &NEVPT2 &END
nr_states = 1

999 continue
key = Get_Ln(LuSpool)
call LeftAd(key)
Line = key
if (Line(1:1) == '*') goto 999
if (Line == ' ') goto 999
call UpCase(Line)
if (Line(1:4) == 'NOPC') goto 1000
if (Line(1:4) == 'STAT') goto 2000
if (Line(1:4) == 'FROZ') goto 3000
if (Line(1:4) == 'SKIP') goto 4000
if (Line(1:4) == 'NOMS') goto 5000
if (Line(1:4) == 'MULT') goto 6000
if (Line(1:4) == 'FILE') goto 7000
if (Line(1:4) == 'RDMR') goto 8000
if (Line(1:4) == 'DIST') goto 9000
if (Line(1:4) == 'END ') Go To 99999
write(u6,*) 'Unidentified key word  : '
call FindErrorLine()
call Quit_OnUserError()

!========= NOPC =============

1000 continue
no_pc = .true.
Go To 999

!========= STAT =============
2000 continue
! read the # of states
key = Get_Ln(LuSpool)
read(key,*,Err=9920) nr_states
allocate(MultGroup%State(nr_states))
do i=1,nr_states
  MultGroup%State(i) = i
end do
Go To 999
!========= FROZ =============
3000 continue
! Read in the information about frozen orbitals
! It can be provided either with a number of frozen orbitals
! or with a list -- in the 2nd case one needs to provide a keyword
! 'Select' immediately after 'Frozen' and provide a number of frozen
! orbitals followed by a list of indexes of the frozen orbitals.
! E.g. either
! Frozen=20 (orbitals from 1-20 are frozen) or
! Frozen
! Select
! 3 1 2 4 -- 3 orbitals 1 2 4 are frozen
key = Get_Ln(LuSpool)
call LeftAd(key)
call UpCase(key)
! If the first line after the keyword contains 'SELECT'
if (key(1:4) == 'SELE') then
  ! read in the number of frozen orbitals, then the list
  frozen_str = Get_Ln(LuSpool)
  call LeftAd(frozen_str)
  read(frozen_str,*,Err=9920) nr_frozen_orb
  if (nr_frozen_orb <= 0) then
    call WarningMessage(2,'Nr of frozen orbitals for selection must be > 0!')
    call Quit_OnUserError()
  end if
  frozen_str = frozen_str(scan(frozen_str,' '):)
  allocate(igelo(nr_frozen_orb)); igelo = 0
  do while (iError < 0)
    read(frozen_str,*,IOStat=iError)(igelo(i),i=1,nr_frozen_orb)
    if (iError > 0) goto 9920
    if (iError < 0) then
      frozen_str = trim(frozen_str)//trim(Get_Ln(LuSpool))
    end if
  end do
else ! read only the number of frozen orbitals and fill the frozen
  ! indexes consecutively from 1 to nr_frozen_orb
  read(key,*,Err=9920) nr_frozen_orb
  if (nr_frozen_orb < 0) then
    call WarningMessage(2,'Nr of frozen orbitals must be >= 0!')
    call Quit_OnUserError()
  else
    if (nr_frozen_orb == 0) then
      write(u6,*) 'Number of frozen orbitals has been set to 0.'
      ! Set it to -1 to signal that frozen orbs have been set forcibly to 0 here
      ! It will be detected in pt2init and reset back to 0
      nr_frozen_orb = -1
    else
      allocate(igelo(nr_frozen_orb))
      do i=1,nr_frozen_orb
        igelo(i) = i
      end do
    end if
  end if
end if
Go To 999
!========= SKIP =============
! Skip the calculation of Koopmans' matrices if requested by the
! SKIP(Koop) keyword
4000 continue
skip_koopro_molcas = .true.
Go To 999
!========= NOMS =============
! Skip the calculation of an effective Hamiltonian (suitable for multistate state-specific calculations)
5000 continue
skip_effective_ham = .true.
Go To 999
!========= MULT =============
! multi-state QD-NEVPT2 calculation requested with the states given below
6000 continue
if (.not. next_non_comment(LuSpool,Line)) goto 9910
read(Line,*) key
call upcase(key)
if (trim(key) == 'ALL') then
  nr_states = 0
else
  read(Line,*,Err=9920,end=9920) nr_states
  if (nr_states <= 0) then
    write(u6,*) ' number of MULT states must be > 0, quitting!'
    call Quit_OnUserError()
  end if
end if
allocate(MultGroup%State(nr_states))
iSplit = scan(Line,' ')
dLine = line(iSplit:)
iError = -1
do while (iError < 0)
  read(dLine,*,IOStat=iError) (MultGroup%State(i),i=1,nr_states)
  if (iError > 0) goto 9920
  if (iError < 0) then
    if (.not. next_non_comment(LuSpool,Line)) goto 9910
    dline = trim(dline)//' '//line
  end if
end do
Go To 999
!========= FILE =============
! Specifiy the name of the reference wfn file for NEVPT2.
7000 continue
if (.not. next_non_comment(LuSpool,Line)) goto 9910
call LeftAd(line)
call fileorb(Line,refwfnfile)
Go To 999
!========= RDMR ============= a.k.a. RDMRead
! Skip calculation of 4-RDM and/or transition 3-RDMs at the beginning of the calculation
! This is useful only for single-file reading, e.g. when your previous calculation has crashed
! and you do not wish to recalculate the full 4-RDM again.
! For distributed RDM calculations see option below
8000 continue
rdm_read = .true.
Go To 999
!========= DIST ============= a.k.a. DistributedRDM
! Skip calculation of RDMs *and* read them from distributed RDM calculation
! Specify a path, where the program will look for subdirectories
! of the format "A-B-C-D", where A,B,C,D are the first four indices of the 4-RDM to be calculated
! Each subdirectory should contain the results of a single calculation in a batch
9000 continue
rdm_distributed = .true.
if (.not. next_non_comment(LuSpool,Line)) goto 9910
read(Line,'(A)') key
rdm_path = trim(key)
Go To 999

! END of Input

9910 continue
call WarningMessage(2,'Premature end of input file.')
goto 9920
9920 continue
call WarningMessage(2,'Read error during input preprocessing.')
call Quit_OnUserError()
call ABEND()

99999 continue

!> make sure the array is allocated for the minimal input
!> &NEVPT2 &END
if (.not. allocated(MultGroup%State)) allocate(MultGroup%State(1)); MultGroup%State(1) = 1

end subroutine rdinput
