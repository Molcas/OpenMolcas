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

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: DIAG(*), EXPLE(*), EXPLV(*)
real(kind=wp), intent(in) :: ONEINT(*), TUVX(*)
integer(kind=iwp), intent(_OUT_) :: ISEL(*)
integer(kind=iwp) :: I, II, IPRLEV, IREOTS, LEXHAM, LOCONE, LW1, LW2, lwscr, MXXSEL, MXXWS, NHEX, NPCNF
real(kind=wp) :: ECORE
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "csfbas.fh"
#include "strnum.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "output_ras.fh"

call Timing(Omega_1,Swatch,Swatch,Swatch)
IPRLEV = IPRLOC(3)

ECORE = Zero
MXXSEL = NSEL
NHEX = NSEL*(NSEL+1)/2

! ALLOCATE LOCAL MEMORY

call GETMEM('IPCNF','ALLO','INTE',LW1,NCNASM(STSYM))
call GETMEM('HONE','ALLO','REAL',LOCONE,NAC**2)
call GETMEM('EXHAM','ALLO','REAL',LEXHAM,NHEX)

! EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE

call TRIEXP(ONEINT,Work(LOCONE),NAC)

! Load the diagonal approximation of the CI Hamiltonian

call Load_H_diag(nConf,DIAG,LuDavid)

! CONSTRUCT THE EXPLICIT HAMILTONIAN

IPRINT = 0
if (IPRLEV == INSANE) IPRINT = 40
call GETMEM('IREOTS','ALLO','INTEGER',IREOTS,NAC)
call GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
call GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)
call GET_IREOTS(IWORK(IREOTS),NAC)
call PHPCSF(Work(LEXHAM),ISEL,iWork(LW1),MXXSEL,Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),STSYM,Work(LOCONE),ECORE,NAC,Work(LW2), &
            NCNASM(STSYM),NAEL+NBEL,NAEL,NBEL,NSEL,NPCNF,DIAG,TUVX,IPRINT,ExFac,IWORK(IREOTS))
if (IPRLEV == INSANE) then
  call Square(Work(LEXHAM),EXPLV,1,NSEL,NSEL)
  call RECPRT('Square Explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)
end if
call GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
call GETMEM('IREOTS','FREE','INTEGER',IREOTS,NAC)
call GETMEM('HONE','FREE','REAL',LOCONE,NAC**2)
call GETMEM('IPCNF','FREE','INTE',LW1,NCNASM(STSYM))

! DIAGONALIZE THE EXPLICIT HAMILTONIAN.

!if (nSel == nConf) then
if (.true.) then
  call DCOPY_(MXXSEL*MXXSEL,[Zero],0,EXPLV,1)
  do I=1,NSEL
    II = I+NSEL*(I-1)
    EXPLV(II) = One
  end do
  !call Jacob(Work(LEXHAM),EXPLV,NSEL,NSEL)
  !# ifdef _DEBUGPRINT_
  !call NIdiag(Work(LEXHAM),EXPLV,NSEL,NSEL)
  !# else
  call NIdiag_new(Work(LEXHAM),EXPLV,NSEL,NSEL)
  !# endif
  call JACORD(Work(LEXHAM),EXPLV,NSEL,NSEL)
  do I=1,NSEL
    EXPLE(I) = Work(LEXHAM-1+I*(I+1)/2)
  end do
else
  call GetMem('ExHscr','Allo','Real',lwscr,nSel)
  call Square(Work(LEXHAM),EXPLV,1,NSEL,NSEL)
  call Eigen_Molcas(NSEL,EXPLV,EXPLE,Work(lwscr))
  call GetMem('ExHscr','Free','Real',lwscr,nSel)
end if
call GETMEM('EXHAM','FREE','REAL',LEXHAM,NHEX)
if (IPRLEV >= INSANE) call IVCPRT('Configurations included in the explicit Hamiltonian',' ',ISEL,NSEL)
if (IPRLEV >= INSANE) call DVCPRT('Eigenvalues of the explicit Hamiltonian',' ',EXPLE,NSEL)
if (IPRLEV >= INSANE) call RECPRT('Eigenvectors of the explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)

call Timing(Omega_2,Swatch,Swatch,Swatch)
Omega_2 = Omega_2-Omega_1
Omega_3 = Omega_3+Omega_2

return

end subroutine EXPLH2
