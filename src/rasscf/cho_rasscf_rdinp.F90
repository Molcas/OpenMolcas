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

subroutine CHO_RASSCF_RDINP(DFonly,LuInput)
!***********************************************************************
!                                                                      *
!  Purpose:   If DFonly, use defaults only.                            *
!             Else, read and process input for Cholesky section        *
!             in RASSCF                                                *
!                                                                      *
!***********************************************************************

use Fock_util_global, only: ALGO, Deco, DensityCheck, dmpK, DoLocK, Estimate, Nscreen, Update
use Cholesky, only: ChFracMem, timings
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp) :: DFonly
integer(kind=iwp) :: LuInput
integer(kind=iwp) :: iChrct, iCLast, Last
character(len=180) :: Key, KWord
real(kind=wp) :: DMPK_DFL
character(len=180), external :: Get_Ln

!**** Algorithms for using Cholesky vectors in RASSCF ******************
!
!   ALGO:
!
!          1  --->  compute fock matrices in AO-basis from vectors
!                   read in reduced sets and transformed on the fly
!                   (PU|VX) integrals are also returned
!
!                   If DoLock=.true. then the Exchange terms are
!                   computed by emplpying the "Local K" scheme
!
!
!          2  --->  compute fock matrices in AO-basis from vectors
!                   read in reduced sets and transformed on the fly
!                   Only the (TU|VX) integrals are returned
!                   while the vectors are directly contracted
!                   with the 2-body density in order to construct
!                   the Q-matrix
!                                                                      *
!***********************************************************************
!                                                                      *
! Default  parameters
#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Zero
#endif

! set some parameters if not specified in ChoInput section
ALGO = 1
DensityCheck = .false.
Deco = .true.
timings = .false.
DoLock = .true.
Nscreen = 10
dmpk = 1.0e-1_wp
Update = .true.
Estimate = .false.
if (DFonly) return

dmpk_dfl = 1.0e-1_wp
!***********************************************************************
!                                                                      *
!-----Process the input

!----------------------------------------------------------------------*
! The big turning point.                                               *
!----------------------------------------------------------------------*
outer: do
  !--------------------------------------------------------------------*
  ! Use Get_Ln to read the lines.                                      *
  !--------------------------------------------------------------------*
  Key = Get_Ln(LuInput)
  Kword = Key
  call UpCase(Kword)
  !--------------------------------------------------------------------*
  ! The keywords and their labels.                                     *
  !--------------------------------------------------------------------*

  if ((KWord(1:1) == '*') .or. (KWord == '')) cycle outer
  select case (KWord(1:4))
    case ('ALGO')
      !                                                                *
      !***** ALGO ******************************************************
      !                                                                *
      !-----Read Cholesky algorithm parameters

      !call Get_F1(1,Eps)
      !call Get_F1(2,rds)
      !call Get_I1(1,ALGO)
      !if (nToken(KWord) > 1) exit outer

      read(LuInput,*) ALGO

      if (ALGO == 1) then
        write(u6,*) 'Default RASSCF algorithm reset to  ',ALGO
        write(u6,*)
      else if (ALGO == 2) then
        write(u6,*) 'Default RASSCF algorithm reset to  ',ALGO
        write(u6,*)
        write(u6,*) ' !!! STILL UNDER DEBUGGING !!! '
      else
        write(u6,*) 'The specified algorithm is not implemented. Option Ignored '
        write(u6,*)
      end if

    case ('LOCK','LK')
      !                                                                *
      !***** LOCK or LK ************************************************
      !                                                                *
      DoLocK = .true.
      !write(u6,*) 'Using Local K scheme for Exchange matrices'

    case ('NOLK')
      !                                                                *
      !***** NOLK ******************************************************
      !                                                                *
      DoLocK = .false.
      !write(u6,*) 'LK screening for Exchange matrices turned off!'

    case ('DMPK')
      !                                                                *
      !***** DMPK ******************************************************
      !                                                                *
      read(LuInput,*) dmpk
      if (dmpk < Zero) then
        write(u6,*) 'OBS! Specified Negative DMPK value. Restore Defaults'
        dmpk = dmpk_dfl
      end if

    case ('NODE')
      !                                                                *
      !***** NODE ******************************************************
      !                                                                *
      Deco = .false.
      write(u6,*) 'Not-Using Cholesky decomposed Inactive density'

    case ('SCRN')
      !                                                                *
      !***** SCRN ******************************************************
      !                                                                *
      read(LuInput,*) Nscreen

    case ('MEMF')
      !                                                                *
      !***** MEMF ******************************************************
      !                                                                *
      read(LuInput,*) ChFracMem

    case ('DCHK')
      !                                                                *
      !***** DCHK ******************************************************
      !                                                                *
      DensityCheck = .true.
      write(u6,*) 'Non-valid option. IGNORED !! '

    case ('TIME')
      !                                                                *
      !***** TIME ******************************************************
      !                                                                *
      timings = .true.

    case ('ESTI')
      !                                                                *
      !***** ESTI ******************************************************
      !                                                                *
      Estimate = .true.
      write(u6,*) 'Diagonal integrals estimated from the current Cholesky vectors'

    case ('UPDA')
      !                                                                *
      !***** UPDA ******************************************************
      !                                                                *
      Update = .true.
      write(u6,*) 'Updating of the true diagonal integrals'

    case ('END','ENDC','ENDO')
      !                                                                *
      !***** END  ******************************************************
      !                                                                *
      !-----End of input
      exit outer

    case default
      !----------------------------------------------------------------*
      ! Control section                                                *
      !----------------------------------------------------------------*
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(u6,*) 'CHO_RASSCF_RDINP Error in keyword.'
      call Quit_OnUserError()

  end select

end do outer
!                                                                      *
!***********************************************************************
!                                                                      *

!write(u6,'(1X,A,I4)') 'Default Cholesky algorithm in RASSCF = ',ALGO
write(u6,*)
if ((ALGO == 2) .and. (DoLocK)) then
  write(u6,*) 'Local K scheme not implemented for the chosen algorithm. LocK keyword ignored !'
  DoLocK = .false.
end if

end subroutine CHO_RASSCF_RDINP
