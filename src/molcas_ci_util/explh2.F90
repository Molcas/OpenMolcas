!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1996, Markus P. Fuelscher                         *
!***********************************************************************

subroutine EXPLH2(DIAG,ONEINT,TUVX,ISEL,EXPLE,EXPLV)
!***********************************************************************
!                                                                      *
!     Compute and diagonalize the explicit Hamiltonian in a            *
!     subspace smaller or identical to nSel. nSel may be chosen        *
!     somewhat smaller in order to avoid selecting only one of         *
!     several degenerate diagonal matrix elements.                     *
!                                                                      *
!     calling arguments:                                               *
!     DIAG    : array of Real                                          *
!               diagonal Hamiltonian                                   *
!     ONEINT  : array of Real                                          *
!               one-electron integrals                                 *
!     TUVX    : array of Real                                          *
!               two-electron integrals                                 *
!     ISEL    : array of integer                                       *
!               index array                                            *
!     EXPLE   : array of Real                                          *
!               eigenvalues of the explicit Hamiltonian                *
!     EXPLV   : array of Real                                          *
!               eigenvectors of the explicit Hamiltonian               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and J. Olsen                                      *
!     University of Lund, Sweden, 1990                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use timers, only: TimeHSel
use rasscf_global, only: NAC, NACPR2
use rasscf_files, only: LUDAVID
use CI_Info, only: CI_St
use rasscf_global, only: NCONF, NSEL
use output_ras, only: IPRLOC
use PrintLevel, only: INSANE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
! the RASSCF CI state (CI_Info)
integer(kind=iwp), parameter :: iState = 1
real(kind=wp), intent(_OUT_) :: DIAG(*), EXPLE(*), EXPLV(*)
real(kind=wp), intent(in) :: ONEINT(*), TUVX(*)
integer(kind=iwp), intent(_OUT_) :: ISEL(*)
integer(kind=iwp) :: I, II, IPRINT, IPRLEV, MXXSEL, NHEX
real(kind=wp) :: dum1, dum2, dum3, ECORE, Time(2)
real(kind=wp), allocatable :: EXHAM(:), Scr(:)

call Timing(Time(1),dum1,dum2,dum3)
IPRLEV = IPRLOC(3)

ECORE = Zero
MXXSEL = NSEL
NHEX = nTri_Elem(NSEL)

! ALLOCATE LOCAL MEMORY

call mma_allocate(EXHAM,NHEX,label='EXHAM')

! Load the diagonal approximation of the CI Hamiltonian

call Load_H_diag(nConf,DIAG,LuDavid)

! CONSTRUCT THE EXPLICIT HAMILTONIAN
! (selection of the subspace and <I|H|J> through the CI interface)

IPRINT = 0
if (IPRLEV == INSANE) IPRINT = 40
call PHPCSF(EXHAM,ISEL,MXXSEL,CI_St(iState)%STSym,ECORE,NSEL,DIAG,ONEINT,nTri_Elem(NAC),TUVX,NACPR2,IPRINT)
if (IPRLEV == INSANE) then
  call Square(EXHAM,EXPLV,1,NSEL,NSEL)
  call RECPRT('Square Explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)
end if

! DIAGONALIZE THE EXPLICIT HAMILTONIAN.

!if (nSel == nConf) then
if (.true.) then
  EXPLV(1:NSEL*NSEL) = Zero
  do I=1,NSEL
    II = I+NSEL*(I-1)
    EXPLV(II) = One
  end do
  !call Jacob(EXHAM,EXPLV,NSEL,NSEL)
  !# ifdef _DEBUGPRINT_
  !call NIdiag(EXHAM,EXPLV,NSEL,NSEL)
  !# else
  call NIdiag_new(EXHAM,EXPLV,NSEL,NSEL)
  !# endif
  call JACORD(EXHAM,EXPLV,NSEL,NSEL)
  do I=1,NSEL
    EXPLE(I) = EXHAM(I*(I+1)/2)
  end do
else
  call mma_allocate(Scr,nSel,label='ExHscr')
  call Square(EXHAM,EXPLV,1,NSEL,NSEL)
  call Eigen_Molcas(NSEL,EXPLV,EXPLE,Scr)
  call mma_deallocate(Scr)
end if
call mma_deallocate(EXHAM)
if (IPRLEV >= INSANE) call IVCPRT('Configurations included in the explicit Hamiltonian',' ',ISEL,NSEL)
if (IPRLEV >= INSANE) call DVCPRT('Eigenvalues of the explicit Hamiltonian',' ',EXPLE,NSEL)
if (IPRLEV >= INSANE) call RECPRT('Eigenvectors of the explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)

call Timing(Time(2),dum1,dum2,dum3)
TimeHSel = TimeHSel+Time(2)-Time(1)

return

end subroutine EXPLH2
