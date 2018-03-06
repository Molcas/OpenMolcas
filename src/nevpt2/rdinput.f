************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*
* MOLCAS wrapper for Celestino Angeli's NEVPT2 code
* rdinput reads input from the MOLCAS input file
*
********************
* This routine is based on the EXPBAS input routine in
* ../src/readinp_expbas.f
* therefore we stick to ugly F77 code and use nice F90 features at the same time! O_o
* One day it should be changed to proper F90 code

      subroutine rdinput(refwfnfile)

      use nevpt2_cfg ! use global variables directly from the NEVPT2 program

      implicit none

      Character(len=*), intent(out) :: refwfnfile

      Character*180  Line, Blank, key, Get_Ln
      Character*9001 dline
      Character*9001 frozen_str

      External Get_Ln, isFreeUnit
      integer LuSpool, isFreeUnit, iError,i, isplit
      logical, external :: next_non_comment

c
c Initial values
c

      refwfnfile = "JOBIPH"
c
      LuSpool=18
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)

      Rewind(LuSpool)
      Call RdNLst(LuSpool,'NEVPT2')
      Blank=' '
      iError=-1

      !> set it to one in order to make sure NEVPT2 runs also with just the basic input
      !> &NEVPT2 &END
      nr_states = 1

  999 Continue
*      Read(LuSpool,'(A)',End=9940) Line
      key =Get_Ln(LuSpool)
      Call LeftAd(key)
      Line = key
      If (Line(1:1).eq.'*' ) Goto 999
      If (Line.eq.Blank ) Goto 999
      Call UpCase(Line)
      If (Line(1:4).eq.'NOPC') Goto 1000
      If (Line(1:4).eq.'STAT') Goto 2000
      If (Line(1:4).eq.'FROZ') Goto 3000
      If (Line(1:4).eq.'SKIP') Goto 4000
      If (Line(1:4).eq.'NOMS') Goto 5000
      If (Line(1:4).eq.'MULT') Goto 6000
      If (Line(1:4).eq.'FILE') Goto 7000
      If (Line(1:4).eq.'END ') Go To 99999
      Write (6,*) 'Unidentified key word  : '
      Call FindErrorLine
      Call Quit_OnUserError()

*========= NOPC =============

 1000 Continue
      no_pc = .true.
      Go To 999

*========= STAT =============
 2000 Continue
      ! read the # of states
      key =Get_Ln(LuSpool)
      Read(key,*,Err=9920) nr_states
      Allocate(MultGroup%State(nr_states))
      do i = 1, nr_states
        MultGroup%State(i) = i
      end do
      Go To 999
*========= FROZ =============
 3000 Continue
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
      key =Get_Ln(LuSpool)
      Call LeftAd(key)
      Call UpCase(key)
      ! If the first line after the keyword contains 'SELECT'
      If(key(1:4).eq.'SELE') then
        ! read in the number of frozen orbitals, then the list
        frozen_str =Get_Ln(LuSpool)
        Call LeftAd(frozen_str)
        read (frozen_str,*,Err=9920) nr_frozen_orb
        if (nr_frozen_orb.le.0) then
           Call WarningMessage(2,'Nr of frozen orbitals for selection'
     &      //' must be > 0!')
           Call Quit_OnUserError
        end if
        frozen_str=frozen_str(scan(frozen_str,' '):)
        allocate(igelo(nr_frozen_orb));igelo = 0
        Do While (iError.lt.0)
          Read(frozen_str,*,IOStat=iError) (igelo(i),i=1,nr_frozen_orb)
          If (iError.gt.0) GoTo 9920
          If (iError.lt.0) Then
            frozen_str=trim(frozen_str)//trim(Get_Ln(LuSpool))
          end if
        End Do
      Else ! read only the number of frozen orbitals and fill the frozen
           ! indexes consecutively from 1 to nr_frozen_orb
        read (key,*,Err=9920) nr_frozen_orb
        if (nr_frozen_orb.lt.0) then
          Call WarningMessage(2,'Nr of frozen orbitals must be >= 0!')
          Call Quit_OnUserError
        else
          if (nr_frozen_orb.eq.0) then
            write (6,*) 'Number of frozen orbitals'
     &       //' has been set to 0.'
            ! Set it to -1 to signal that frozen orbs have been set forcibly to 0 here
            ! It will be detected in pt2init and reset back to 0
            nr_frozen_orb = -1
          else
            allocate(igelo(nr_frozen_orb))
            do i=1,nr_frozen_orb
              igelo(i)=i
            end do
          end if
        end if
      end if
      Go To 999
*========= SKIP =============
      ! Skip the calculation of Koopmans' matrices if requested by the
      ! SKIP(Koop) keyword
 4000 Continue
      skip_koopro_molcas = .true.
      Go To 999
*========= NOMS =============
      ! Skip the calculation of an effective Hamiltonian (suitable for multistate state-specific calculations)
 5000 Continue
      skip_effective_ham = .true.
      Go To 999
*========= MULT =============
      ! multi-state QD-NEVPT2 calculation requested with the states given below
 6000 Continue
      If(.NOT.next_non_comment(LuSpool,Line)) GoTo 9910
      Read(Line,*) key
      call upcase(key)
      If (trim(key)=='ALL') Then
        nr_states = 0
      Else
        Read(Line,*,Err=9920,End=9920) nr_states
        If (nr_states.le.0) Then
          Write(6,*)' number of MULT states must be > 0, quitting!'
          Call Quit_OnUserError
        End If
      End If
      Allocate(MultGroup%State(nr_states))
      iSplit = scan(Line,' ')
      dLine = line(iSplit:)
      iError = -1
      Do While (iError.lt.0)
        Read(dLine,*,IOStat=iError)
     &    (MultGroup%State(i), i=1,nr_states)
        If (iError.gt.0) GoTo 9920
        If (iError.lt.0) Then
          If(.NOT.next_non_comment(LuSpool,Line)) GoTo 9910
           dline = trim(dline) // ' ' // line
        End If
      End Do
      Go To 999
*========= FILE =============
      ! Specifiy the name of the reference wfn file for NEVPT2.
 7000 Continue
      If(.NOT.next_non_comment(LuSpool,Line)) GoTo 9910
      Call LeftAd(line)
      Read(Line,*,Err=9920,End=9920) refwfnfile
      Go To 999
c
c END of Input
c

9910  Continue
      Call WarningMessage(2,'Premature end of input file.')
      GoTo 9920
9920  CONTINUE
      Call WarningMessage(2,'Read error during input preprocessing.')
      Call Quit_OnUserError()
      CALL ABEND()

99999 Continue

      !> make sure the array is allocated for the minimal input
      !> &NEVPT2 &END
      if(.not.allocated(MultGroup%State))
     &  allocate(MultGroup%State(1)); MultGroup%State(1) = 1

      end subroutine rdinput
