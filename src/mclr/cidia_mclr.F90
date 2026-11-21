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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine CIDIA_MCLR(iSym,ralp)

use ipPage, only: ipclose, ipget, ipin, ipnout, W
use CI_Interfaces_MCLR, only: CI_MCLR_Diag
use MCLR_Data, only: FIMO, INT2, ipCI, ipDia, NCSASM, nexp, nexp_max, nGP, XISPSM
use input_mclr, only: ERASSCF, nCSF, PotNuc, rIn_Ene, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(out) :: ralp
integer(kind=iwp) :: iC, ipDIAI, nD, nP1, nP2, nQ, nSD
real(kind=wp) :: ECAS
real(kind=wp), allocatable :: Q(:)
real(kind=wp), external :: DDot_

! This is just a interface to hide Jeppe from the rest of the world
! we dont want to let world see the work of the Danish
! (I hope he never reads that)
! Anyway concerning the CSF/SD stuff.
! If we work with spin dependent perturbations
! we never use CSF's (to complicated), instead we use
! SD in all parts of the program,
! otherwise we will switch to SD representation in this routine

! The CSF diagonal through the CI interface (CI_MCLR_Diag)
nsd = max(ncsf(isym),nint(XISPSM(ISYM,1)))
ipdiai = ipget(nsd)
call ipin(ipdiai)
call CI_MCLR_Diag(iState,iSym,size(FIMO),FIMO,size(INT2),INT2,rin_ene+potnuc,ncsf(iSym),W(ipdiai)%A(1:ncsf(iSym)))
nD = NCSASM(ISYM)

! Calculate explicit part of hamiltonian

np2 = min(nd,nexp_max)
np1 = 0
nq = 0
if (np2 /= 0) then
  call ipnout(ipdiai)
  call ipin(ipdiai)
  call h0(W(ipdiai)%A,np1,nexp_max,nq,isym,nexp,TimeDep)
else
  nexp = 0
end if

ECAS = ERASSCF(1)
call ipin(ipdiai)
do iC=1,nD
  if ((W(ipdiai)%A(ic)-ECAS) /= Zero) then
    W(ipdiai)%A(iC) = One/(W(ipdiai)%A(iC)-ECAS)
  else
    W(ipdiai)%A(iC) = 1.0e5_wp
  end if
end do
!             -1
! ralp=<0|(H-E) |0>

call mma_allocate(Q,nD,Label='Q')
Q(:) = Zero

call ipin(ipCI)
call ExpHinvv(W(ipdiai)%A,W(ipCI)%A,Q,Zero,One)

ralp = DDOT_(nD,W(ipCI)%A,1,Q,1)
if (NGP) then
  call MKP1INV(W(ipdiai)%A)
  call MKCIPRE()
end if
call ipnout(ipdiai)
call mma_deallocate(Q)

ipdia = ipdiai

end subroutine CIDIA_MCLR
