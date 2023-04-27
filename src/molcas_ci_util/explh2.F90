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
!     DIAG    : array of Real*8                                        *
!               diagonal Hamiltonian                                   *
!     ONEINT  : array of Real*8                                        *
!               one-electron integrals                                 *
!     TUVX    : array of Real*8                                        *
!               two-electron integrals                                 *
!     ISEL    : array of integer                                       *
!               index array                                            *
!     EXPLE   : array of Real*8                                        *
!               eigenvalues of the explicit Hamiltonian                *
!     EXPLV   : array of Real*8                                        *
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

use csfbas, only: CONF, KDFTP, KDTOC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: DIAG(*), EXPLE(*), EXPLV(*)
real(kind=wp), intent(in) :: ONEINT(*), TUVX(*)
integer(kind=iwp), intent(_OUT_) :: ISEL(*)
integer(kind=iwp) :: I, II, IPRLEV, MXXSEL, MXXWS, NHEX, NPCNF
real(kind=wp) :: dum1, dum2, dum3, ECORE
integer(kind=iwp), allocatable :: CNF(:), IREOTS(:)
real(kind=wp), allocatable :: EXHAM(:), HONE(:,:), Scr(:)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "strnum.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "output_ras.fh"

call Timing(Omega_1,dum1,dum2,dum3)
IPRLEV = IPRLOC(3)

ECORE = Zero
MXXSEL = NSEL
NHEX = NSEL*(NSEL+1)/2

! ALLOCATE LOCAL MEMORY

call mma_allocate(CNF,NCNASM(STSYM),label='IPCNF')
call mma_allocate(HONE,NAC,NAC,label='HONE')
call mma_allocate(EXHAM,NHEX,label='EXHAM')

! EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE

call SQUARE(ONEINT,HONE,NAC,1,NAC)

! Load the diagonal approximation of the CI Hamiltonian

call Load_H_diag(nConf,DIAG,LuDavid)

! CONSTRUCT THE EXPLICIT HAMILTONIAN

IPRINT = 0
if (IPRLEV == INSANE) IPRINT = 40
call mma_allocate(IREOTS,NAC,label='IREOTS')
call mma_maxDBLE(MXXWS)
call mma_allocate(Scr,MXXWS,label='EXHSCR')
call GET_IREOTS(IREOTS,NAC)
call PHPCSF(EXHAM,ISEL,CNF,MXXSEL,Work(KDTOC),iWork(KDFTP),CONF,STSYM,HONE,ECORE,NAC,Scr,NCNASM(STSYM),NAEL+NBEL,NAEL,NBEL,NSEL, &
            NPCNF,DIAG,TUVX,IPRINT,ExFac,IREOTS)
if (IPRLEV == INSANE) then
  call Square(EXHAM,EXPLV,1,NSEL,NSEL)
  call RECPRT('Square Explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)
end if
call mma_deallocate(IREOTS)
call mma_deallocate(Scr)
call mma_deallocate(CNF)
call mma_deallocate(HONE)

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

call Timing(Omega_2,dum1,dum2,dum3)
Omega_2 = Omega_2-Omega_1
Omega_3 = Omega_3+Omega_2

return

end subroutine EXPLH2
