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
use text_file, only: extend_line, next_non_comment
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(out) :: refwfnfile
integer(kind=iwp) :: LuSpool, iError, i, isplit
character(len=180) :: Line, key
character(len=9001) :: frozen_str
character(len=:), allocatable :: dLine, Line2
integer(kind=iwp), external :: isFreeUnit
character(len=180), external :: Get_Ln

! Initial values

refwfnfile = 'JOBIPH'

LuSpool = isFreeUnit(18)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'NEVPT2')

!> set it to one in order to make sure NEVPT2 runs also with just the basic input
!> &NEVPT2 &END
nr_states = 1

do
  key = Get_Ln(LuSpool)
  Line = adjustl(key)
  if (Line(1:1) == '*') cycle
  if (Line == ' ') cycle
  call UpCase(Line)
  select case (Line(1:4))
    case ('NOPC')
      !========= NOPC =============
      no_pc = .true.

    case ('STAT')
      !========= STAT =============
      ! read the # of states
      key = Get_Ln(LuSpool)
      read(key,*,iostat=iError) nr_states
      if (iError > 0) call error(0)
      ! using standard allocate and deallocate because MultGroup%State
      ! is deallocated somewhere in the external library
      if (allocated(MultGroup%State)) deallocate(MultGroup%State)
      allocate(MultGroup%State(nr_states))
      do i=1,nr_states
        MultGroup%State(i) = i
      end do

    case ('FROZ')
      !========= FROZ =============
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
      call UpCase(key)
      key = adjustl(key)
      ! If the first line after the keyword contains 'SELECT'
      if (key(1:4) == 'SELE') then
        ! read in the number of frozen orbitals, then the list
        frozen_str = Get_Ln(LuSpool)
        read(frozen_str,*,iostat=iError) nr_frozen_orb
        if (iError > 0) call error(0)
        if (nr_frozen_orb <= 0) then
          call WarningMessage(2,'Nr of frozen orbitals for selection must be > 0!')
          call Quit_OnUserError()
        end if
        frozen_str = frozen_str(scan(frozen_str,' '):)
        ! using standard allocate because igelo
        ! is deallocated somewhere in the external library
        allocate(igelo(nr_frozen_orb))
        igelo(:) = 0
        iError = -1
        do while (iError < 0)
          read(frozen_str,*,iostat=iError) (igelo(i),i=1,nr_frozen_orb)
          if (iError > 0) call error(0)
          if (iError < 0) then
            frozen_str = trim(frozen_str)//trim(Get_Ln(LuSpool))
          end if
        end do
      else ! read only the number of frozen orbitals and fill the frozen
        ! indexes consecutively from 1 to nr_frozen_orb
        read(key,*,iostat=iError) nr_frozen_orb
        if (iError > 0) call error(0)
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
            ! using standard allocate because igelo
            ! is deallocated somewhere in the external library
            allocate(igelo(nr_frozen_orb))
            do i=1,nr_frozen_orb
              igelo(i) = i
            end do
          end if
        end if
      end if

    case ('SKIP')
      !========= SKIP =============
      ! Skip the calculation of Koopmans' matrices if requested by the
      ! SKIP(Koop) keyword
      skip_koopro_molcas = .true.

    case ('NOMS')
      !========= NOMS =============
      ! Skip the calculation of an effective Hamiltonian (suitable for multistate state-specific calculations)
      skip_effective_ham = .true.

    case ('MULT')
      !========= MULT =============
      ! multi-state QD-NEVPT2 calculation requested with the states given below
      if (.not. next_non_comment(LuSpool,Line2)) call error(1)
      read(Line2,*) key
      call upcase(key)
      if (trim(key) == 'ALL') then
        nr_states = 0
      else
        read(Line2,*,iostat=iError) nr_states
        if (iError /= 0) call error(0)
        if (nr_states <= 0) then
          write(u6,*) ' number of MULT states must be > 0, quitting!'
          call Quit_OnUserError()
        end if
      end if
      ! using standard allocate and deallocate because MultGroup%State
      ! is deallocated somewhere in the external library
      if (allocated(MultGroup%State)) deallocate(MultGroup%State)
      allocate(MultGroup%State(nr_states))
      iSplit = scan(Line2,' ')
      call mma_allocate (dLine,len(Line),label='dLine')
      dLine = line2(iSplit:)
      iError = -1
      do while (iError < 0)
        read(dLine,*,iostat=iError) (MultGroup%State(i),i=1,nr_states)
        if (iError > 0) call error(0)
        if (iError < 0) then
          if (.not. next_non_comment(LuSpool,Line2)) call error(1)
          call extend_line(dLine,Line)
        end if
      end do
      call mma_deallocate (dLine)

    case ('FILE')
      !========= FILE =============
      ! Specifiy the name of the reference wfn file for NEVPT2.
      if (.not. next_non_comment(LuSpool,Line2)) call error(1)
      line2(:) = adjustl(line2)
      call fileorb(Line2,refwfnfile)

    case ('RDMR')
      !========= RDMR ============= a.k.a. RDMRead
      ! Skip calculation of 4-RDM and/or transition 3-RDMs at the beginning of the calculation
      ! This is useful only for single-file reading, e.g. when your previous calculation has crashed
      ! and you do not wish to recalculate the full 4-RDM again.
      ! For distributed RDM calculations see option below
      rdm_read = .true.

    case ('DIST')
      !========= DIST ============= a.k.a. DistributedRDM
      ! Skip calculation of RDMs *and* read them from distributed RDM calculation
      ! Specify a path, where the program will look for subdirectories
      ! of the format "A-B-C-D", where A,B,C,D are the first four indices of the 4-RDM to be calculated
      ! Each subdirectory should contain the results of a single calculation in a batch
      rdm_distributed = .true.
      if (.not. next_non_comment(LuSpool,Line2)) call error(1)
      read(Line2,'(A)') key
      rdm_path = trim(key)

    case ('END ')
      exit

    case default
      write(u6,*) 'Unidentified key word  : '
      call FindErrorLine()
      call Quit_OnUserError()

  end select

end do
! END of Input

if (allocated(Line2)) call mma_deallocate(Line2)

!> make sure the array is allocated for the minimal input
!> &NEVPT2 &END
if (.not. allocated(MultGroup%State)) then
  ! using standard allocate because MultGroup%State
  ! is deallocated somewhere in the external library
  allocate(MultGroup%State(1))
  MultGroup%State(:) = 1
end if

contains

subroutine error(code)

  integer(kind=iwp), intent(in) :: code

  if (code == 1) call WarningMessage(2,'Premature end of input file.')
  call WarningMessage(2,'Read error during input preprocessing.')
  call Quit_OnUserError()
  call ABEND()

end subroutine error

end subroutine rdinput
