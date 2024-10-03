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

!#define _DEBUGPRINT_
subroutine CHO_SCF_RDINP(DFonly,LuSpool)
!***********************************************************************
!
!  Purpose:   If DFonly, use defaults only.
!             Else, read and process input for Cholesky section in SCF
!
!***********************************************************************

use Fock_util_global, only: Deco, DensityCheck, Estimate, Update
use Cholesky, only: ChFracMem, timings
use InfSCF, only: ALGO, dmpk, nScreen, ReOrd
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: DFonly
integer(kind=iwp), intent(in) :: LuSpool
#include "print.fh"
integer(kind=iwp) :: i, iChrct, iPrint, iRout, jRout, Last, n
real(kind=wp) :: dmpk_dfl
character(len=180) :: Key, KWord
character(len=*), parameter :: SECNAM = 'CHO_SCF_RDINP'
integer(kind=iwp), external :: iCLast
character(len=180), external :: Get_Ln

iRout = 1
iPrint = nPrint(iRout)
!                                                                      *
!**** Algorithms for using Cholesky vectors in SCF *********************
!                                                                      *
!   ALGO:
!          0  --->  Integrals are regenerated on the fly
!                   from a set of Cholesky vectors resorted on disk
!                   (Used only for debugging! Amazingly slow!)
!
!          1  --->  The resorted Cholesky vectors are used directly
!                   by the Fock matrix builder routines and contracted
!                   with the proper density matrices. Uses
!                   vectors resorted either on disk or on the fly
!                   (Used only for debugging! Very slow!)
!
!          2  --->  As in option 1 but using the MO-basis transformed
!                   vectors for computing the exchange term
!
!          3  --->  As in option 2 but using the MO-basis vectors
!                   transformed directly in reduced sets
!
!          4  --->  Local-exchange (LK) algorithm for the exchange term
!                                                                      *
!***********************************************************************
!                                                                      *
! Default  parameters

#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Half
#endif

! set some parameters if not specified in ChoInput section
ALGO = 4
REORD = .false.
DECO = .true.
DensityCheck = .false.
timings = .false.
NSCREEN = 10    ! default screening interval (# of red sets)
dmpk = 0.045_wp ! default damping of the screening threshold
Estimate = .false.
Update = .true.
if (DFonly) return

dmpk_dfl = dmpk
!                                                                      *
!***********************************************************************
!                                                                      *
iPrint = 5
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the input

!----------------------------------------------------------------------*
! The big turning point.                                               *
!----------------------------------------------------------------------*
do
  !--------------------------------------------------------------------*
  ! Use Get_Ln to read the lines.                                      *
  !--------------------------------------------------------------------*
  Key = Get_Ln(LuSpool)
  Kword = Key
  call UpCase(Kword)
  !--------------------------------------------------------------------*
  ! The keywords and their labels.                                     *
  !--------------------------------------------------------------------*

  if ((KWord(1:1) == '*') .or. (KWord == '')) cycle
  select case (KWord(1:4))
    case ('ALGO')
      !                                                                *
      !***** ALGO ******************************************************
      !                                                                *
      ! Read Cholesky algorithm parameters

      read(LuSpool,*) ALGO

#     ifdef _DEBUGPRINT_
      select case (ALGO)
        case (0)
          write(u6,*) 'Integral regeneration from Cholesky vectors reordered on disk'
        case (1)
          write(u6,*) 'Density-based Cholesky. Default reorder: on the fly'
        case (2)
          write(u6,*) 'MO-based-Exchange Cholesky. Default reorder: on the fly'
        case (3)
          write(u6,*) 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
        case (4)
          write(u6,*) 'Local-Exchange (LK) algorithm.'
      end select
      write(u6,*)
#     endif

    case ('REOR')
      !                                                                *
      !***** REOR ******************************************************
      !                                                                *
      REORD = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Vectors reordered on DISK'
      write(u6,*)
#     endif

    case ('NODE')
      !                                                                *
      !***** NODE ******************************************************
      !                                                                *
      DECO = .false.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Not-Using Decomposed density matrix'
      write(u6,*)
#     endif

    case ('DCHK')
      !                                                                *
      !***** DCHK ******************************************************
      !                                                                *
      DensityCheck = .true.

    case ('TIME')
      !                                                                *
      !***** TIME ******************************************************
      !                                                                *
      timings = .true.

    case ('SCRN')
      !                                                                *
      !***** SCRN ******************************************************
      !                                                                *
      read(LuSpool,*) NSCREEN

    case ('DMPK')
      !                                                                *
      !***** DMPK ******************************************************
      !                                                                *
      read(LuSpool,*) dmpk
      if (dmpk < Zero) then
        write(u6,*) 'OBS! Specified Negative DMPK value. Restore Defaults'
        dmpk = dmpk_dfl
      end if

    case ('UPDA')
      !                                                                *
      !***** UPDA ******************************************************
      !                                                                *
      Update = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Local-K with updating of the true diagonals'
      write(u6,*)
#     endif

    case ('ESTI')
      !                                                                *
      !***** ESTI ******************************************************
      !                                                                *
      Estimate = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Local-K with evaluation of the diagonals from the current vec'
      write(u6,*)
#     endif

    case ('LOCK','LK  ')
      !                                                                *
      !***** LOCK or LK ************************************************
      !                                                                *
      algo = 4
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Local-Exchange (LK) algorithm.'
      write(u6,*)
#     endif

    case ('NOLK')
      !                                                                *
      !***** NoLK ******************************************************
      !                                                                *
      algo = 3
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Local-Exchange (LK) screening turned off! '
      write(u6,*)
#     endif

    case ('MEMF')
      !                                                                *
      !***** MemF ******************************************************
      !                                                                *
      read(LuSpool,*) ChFracMem

    case ('PRIN')
      !                                                                *
      !***** PRIN ******************************************************
      !                                                                *
      ! Print level

      Key = Get_Ln(LuSpool)
      KWord = Key
      call Get_I1(1,n)
      do i=1,n
        KWord = Get_Ln(LuSpool)
        call Get_I1(1,jRout)
        call Get_I1(2,iPrint)
        nPrint(jRout) = iPrint
      end do

    case ('ENDC','END ','ENDO')
      !                                                                *
      !***** ENDOFchoinput  ********************************************
      !                                                                *
      exit

    case default
      !----------------------------------------------------------------*
      ! Control section                                                *
      !----------------------------------------------------------------*
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(u6,*) SECNAM,' Error in keyword.'
      call Quit_OnUserError()
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine CHO_SCF_RDINP
