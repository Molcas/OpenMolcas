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

subroutine CHO_RASSI_RDINP(DFonly,LuSpool)
!***********************************************************************
!
!  Purpose:   If DFonly, use defaults only.
!             Else, read and process input for Cholesky section
!             in RASSI
!
!***********************************************************************

use Fock_util_global, only: Deco, Estimate, PseudoChoMOs, Update
use Cholesky, only: timings
use cntrl, only: ALGO, Nscreen, dmpk
use rassi_data, only: CHFRACMEM
use Constants, only: Zero
use Definitions, only: wp, u6

implicit none
logical DFonly
integer LuSpool
character(len=180) KWord, Key
character(len=180), external :: Get_Ln
integer iChrct, last
integer, external :: iCLast
real*8 dmpk_dfl

!**** Algorithms for using Cholesky vectors in RASSI ******************
!
!   ALGO:
!
!          1  --->  compute fock matrices in AO-basis from vectors
!                   read in reduced sets and transformed on the fly.
!                   (TU|VX) integrals are also returned
!
!          2  --->  compute fock matrices in AO-basis from vectors
!                   read in reduced sets and transformed on the fly
!                   (TU|VX) integrals are also returned.
!                   Exchange contributions are computed using the
!                   "Local K" scheme
!                                                                      *
!***********************************************************************
!                                                                      *
! Default  parameters
#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Zero
#endif
ALGO = 2
timings = .false.
dmpk = 0.1_wp
Nscreen = 10  ! to be reset to 0 if doing RI calculation
Deco = .true.
Update = .true.
Estimate = .false.

if (.not. DFonly) then

  ! set some parameters if not specified in ChoInput section
  PseudoChoMOs = .false.
  dmpk_dfl = 0.1_wp
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Process the input

  !--------------------------------------------------------------------*
  do
    !------------------------------------------------------------------*
    ! Use Get_Ln to read the lines.                                    *
    !------------------------------------------------------------------*
    Key = Get_Ln(LuSpool)
    Kword = Key
    call UpCase(Kword)
    !------------------------------------------------------------------*
    ! The keywords and their labels.                                   *
    !------------------------------------------------------------------*
    if (KWord(1:1) == '*') cycle
    if (KWord == '') cycle
    select case (Kword(1:4))

      case ('ALGO')
        !                                                              *
        !***** ALGO ****************************************************
        !                                                              *
        ! Read Cholesky algorithm parameters

        ! call Get_F1(1,Eps)
        ! call Get_F1(2,rds)
        ! call Get_I1(1,ALGO)
        ! if (nToken(KWord) > 1) call abend()

        read(LuSpool,*) ALGO
        if (ALGO == 1) then
          write(u6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
          write(u6,*)
        else if (ALGO == 2) then
          write(u6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
          write(u6,*)
        else
          write(u6,*) 'The specified algorithm is not implemented.'
          write(u6,*) 'Option Ignored'
          write(u6,*)
        end if

      case ('LOCK')
        algo = 2
        write(u6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        write(u6,*) 'Using Local K scheme for Exchange matrices'

      case ('NOLK')
        algo = 1
        write(u6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        write(u6,*) 'Local K scheme for Exchange matrices turned off !'

      case ('SCRN')
        read(LuSpool,*) Nscreen

      case ('DMPK')
        read(LuSpool,*) DMPK
        if (dmpk < Zero) then
          write(u6,*) 'OBS! Specified Negative DMPK value. Restore Defaults'
          dmpk = dmpk_dfl
        end if

      case ('TIME')
        timings = .true.

      case ('UPDA')
        Update = .true.
        write(u6,*) 'Local-K with updating of the true diagonals'
        write(u6,*)

      case ('ESTI')
        Estimate = .true.
        write(u6,*) 'Local-K with evaluation of the diagonals from the current vec'
        write(u6,*)

      case ('MEMF')
        read(LuSpool,*) ChFracMem

      case ('NODE')
        if (PseudoChoMOs) then
          write(u6,*) ' The keyword NODEcompose is incompatible with the previously specified keyword PSEUdo. '// &
                      'NODEcompose will be ignored'
        else
          Deco = .false.
          write(u6,*) 'Canonical inactive orbitals used in LK CD-RASSI.'
        end if
        write(u6,*)

      case ('PSEU')
        if (.not. Deco) then
          write(u6,*) ' The keyword PSEUdo is incompatible with the previously specified keyword NODEcompose. PSEUdo will be '// &
                      'ignored'
        else
          PseudoChoMOs = .true.
          write(u6,*) 'Pseudo Cholesky orbitals used in LK CD-RASSI.'
        end if
        write(u6,*)

      case ('ENDC','END ','ENDO')
        exit

      case default
        iChrct = len(KWord)
        Last = iCLast(KWord,iChrct)
        write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
        write(u6,*) 'CHO_RASSI_RDINP: Error in keyword.'
        call AbEnd()

    end select

  end do

end if

end subroutine CHO_RASSI_RDINP
