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

subroutine CHO_CASPT2_RDINP(DFonly,LuSpool)
!***********************************************************************
!
!  Purpose:   If DFonly, use defaults only.
!             Else, read and process input for Cholesky section
!             in CASPt2
!
!***********************************************************************

use Fock_util_global, only: ALGO, Deco, DensityCheck, REORD
use definitions, only: iwp, u6
use Cholesky, only: timings
use ChoCASPT2, only: iAlGO

implicit none
logical(kind=iwp), intent(in) :: DFonly
integer(kind=iwp), intent(in) :: LuSpool
character(Len=180) KWord, Key
character(Len=180), external :: Get_Ln
character(len=16), parameter :: SECNAM = 'CHO_CASPT2_RDINP'
integer(kind=iwp) iChrct, Last
integer(kind=iwp), external :: iCLast

!                                                                      *
!***********************************************************************
!                                                                      *
! Algorithms for generating MO integrals in CASPT2
!
!    iALGO :
!           0  --> MOLINT file is generated from the
!                  transformed Cholesky vectors. Both "Coulomb" and
!                  "Exchange(1,2)" integrals are computed and stored
!                  on disk
!! Is iALGO = 0 really working?
!
!           1  --> Only the "Exchange" integrals are computed and
!                  combined directly in order to compute the RHS
!                  of the caspt2 equations. The latter is then
!                  stored on disk. The AO Fock matrix is
!                  computed during the MO transformation of the
!                  vectors and it is stored on disk
!
!***********************************************************************
!
!**** Algorithms for using Cholesky vectors in Fock matrix generation **
!
!   ALGO:
!      0  --->  Integrals are regenerated on the fly
!               from a set of Cholesky vectors resorted on disk
!
!      1  --->  The resorted Cholesky vectors are used directly
!               by the Fock matrix builder routines and contracted
!               with the proper density matrices. Uses
!               vectors resorted either on disk or on the fly
!
!      2  --->  As in option 1 but using the MO-basis transformed
!               vectors for computing the exchange term
!
!                                                                      *
!***********************************************************************
!                                                                      *
! Default  parameters

if (DFonly) then
  iAlGO = 1
  ALGO = 2
  REORD = .false.
  DECO = .true.
  DensityCheck = .false.
  timings = .false.
  return
end if

! set some parameters if not specified in ChoInput section
iAlGO = 1
ALGO = 2
REORD = .false.
DECO = .true.
DensityCheck = .false.
timings = .false.
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

  if (KWord(1:1) == '*') cycle
  if (KWord == '') cycle

  select case (KWord(1:4))

    case ('ALGO')
      ! Read Cholesky algorithm parameters
      read(LuSpool,*) ALGO

    case ('IALG')
      read(LuSpool,*) iALGO

    case ('REOR')
      REORD = .true.
      write(u6,*) 'Vectors reordered on FILE'
      write(u6,*)

    case ('DECO')
      DECO = .true.
      write(u6,*) 'Decomposed densty matrix'
      write(u6,*)

    case ('TIME')
      timings = .true.

    case ('DCHK')
      DensityCheck = .true.

    case ('END ','ENDO')
      ! End of input
      return
    case Default
      !----------------------------------------------------------------*
      ! Control section                                                *
      !----------------------------------------------------------------*
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(u6,*) SECNAM,' Error in keyword.'
      call ABEND()
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine CHO_CASPT2_RDINP
