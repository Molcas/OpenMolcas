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
use output_ras, only: LF

implicit none
logical DFonly
integer LuInput
character(len=180) KWord, Key
character(len=180), external :: Get_Ln
character(len=16), parameter :: SECNAM = 'CHO_RASSCF_RDINP'
integer iChrct, Last, iCLast
real*8 DMPK_DFL

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
ChFracMem = 0.3d0
#else
ChFracMem = 0.0d0
#endif

! set some parameters if not specified in ChoInput section
ALGO = 1
DensityCheck = .false.
Deco = .true.
timings = .false.
DoLock = .true.
Nscreen = 10
dmpk = 1.0d-1
Update = .true.
Estimate = .false.
if (DFonly) goto 999  !return flag

dmpk_dfl = 1.0d-1
!***********************************************************************
!                                                                      *
!-----Process the input

!----------------------------------------------------------------------*
! The big turning point.                                               *
!----------------------------------------------------------------------*
1000 continue
!----------------------------------------------------------------------*
! Use Get_Ln to read the lines.                                        *
!----------------------------------------------------------------------*
Key = Get_Ln(LuInput)
Kword = Key
call UpCase(Kword)
!----------------------------------------------------------------------*
! The keywords and their labels.                                       *
!----------------------------------------------------------------------*

if (KWord(1:1) == '*') Go To 1000
if (KWord == '') Go To 1000
if (KWord(1:4) == 'ALGO') Go To 900
if (KWord(1:4) == 'LOCK') Go To 910
if (KWord(1:4) == 'LK  ') Go To 910
if (KWord(1:4) == 'NOLK') Go To 915
if (KWord(1:4) == 'DMPK') Go To 920
if (KWord(1:4) == 'NODE') Go To 930
if (KWord(1:4) == 'SCRN') Go To 940
if (KWord(1:4) == 'MEMF') Go To 950
if (KWord(1:4) == 'DCHK') Go To 820
if (KWord(1:4) == 'TIME') Go To 830
if (KWord(1:4) == 'ESTI') Go To 840
if (KWord(1:4) == 'UPDA') Go To 850
if (KWord(1:4) == 'ENDC') Go To 998
if (KWord(1:4) == 'END ') Go To 998
if (KWord(1:4) == 'ENDO') Go To 998

!----------------------------------------------------------------------*
! Control section                                                      *
!----------------------------------------------------------------------*
iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(LF,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
write(LF,*) SECNAM,' Error in keyword.'
call Quit_OnUserError()
!                                                                      *
!***** ALGO ************************************************************
!                                                                      *
!-----Read Cholesky algorithm parameters

900 continue
!call Get_F1(1,Eps)
!call Get_F1(2,rds)
!call Get_I1(1,ALGO)
!if (nToken(KWord) > 1) goto 988

read(LuInput,*) ALGO

if (ALGO == 1) then
  write(LF,*) 'Default RASSCF algorithm reset to  ',ALGO
  write(LF,*)
else if (ALGO == 2) then
  write(LF,*) 'Default RASSCF algorithm reset to  ',ALGO
  write(LF,*)
  write(LF,*) ' !!! STILL UNDER DEBUGGING !!! '
else
  write(LF,*) 'The specified algorithm is not implemented. Option Ignored '
  write(LF,*)
end if

Go To 1000
!                                                                      *
!***** LOCK or LK ******************************************************
!                                                                      *
910 continue
DoLocK = .true.
!write(LF,*) 'Using Local K scheme for Exchange matrices'

Go To 1000
!                                                                      *
!***** NOLK ************************************************************
!                                                                      *
915 continue
DoLocK = .false.
!write(LF,*) 'LK screening for Exchange matrices turned off!'
!
Go To 1000
!                                                                      *
!***** DMPK ************************************************************
!                                                                      *
920 continue
read(LuInput,*) dmpk
if (dmpk < 0.0d0) then
  write(6,*) 'OBS! Specified Negative DMPK value. Restore Defaults'
  dmpk = dmpk_dfl
end if

Go To 1000
!                                                                      *
!***** NODE ************************************************************
!                                                                      *
930 continue
Deco = .false.
write(LF,*) 'Not-Using Cholesky decomposed Inactive density'

Go To 1000
!                                                                      *
!***** SCRN ************************************************************
!                                                                      *
940 continue
read(LuInput,*) Nscreen

Go To 1000
!                                                                      *
!***** MEMF ************************************************************
!                                                                      *
950 continue
read(LuInput,*) ChFracMem

Go To 1000
!                                                                      *
!***** DCHK ************************************************************
!                                                                      *
820 continue
DensityCheck = .true.
write(LF,*) 'Non-valid option. IGNORED !! '

Go To 1000
!                                                                      *
!***** ESTI ************************************************************
!                                                                      *
840 continue
Estimate = .true.
write(LF,*) 'Diagonal integrals estimated from the current Cholesky vectors'

Go To 1000
!                                                                      *
!***** UPDA ************************************************************
!                                                                      *
850 continue
Update = .true.
write(LF,*) 'Updating of the true diagonal integrals'

Go To 1000
!                                                                      *
!***** TIME ************************************************************
!                                                                      *
830 continue
timings = .true.

Go To 1000
!                                                                      *
!***** END  ************************************************************
!                                                                      *
!-----End of input

998 continue
!                                                                      *
!***********************************************************************
!                                                                      *
999 continue
!write(LF,'(1X,A,I4)') 'Default Cholesky algorithm in RASSCF = ',ALGO
write(LF,*)
if (ALGO == 2) then
  write(LF,*) 'Local K scheme not implemented for the chosen algorithm. LocK keyword ignored !'
  DoLocK = .false.
end if

end subroutine CHO_RASSCF_RDINP
