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

use Str_Info, only: CNSM
use ipPage, only: ipclose, ipget, ipin, ipnout, W
use MCLR_Data, only: ICISTR, ipCI, ipDia, NCNATS, NCPCNT, NCSASM, NDPCNT, nexp, nexp_max, nGP, NOCSF, NTYP, XISPSM
use input_mclr, only: ERASSCF, nCSF, PotNuc, rIn_Ene, State_Sym, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(out) :: ralp
integer(kind=iwp) :: i, iAMCmp, iC, iPDCSFI, ipDIAI, iPDSDI, iSM(1), iSPC(1), nD, nP1, nP2, nQ, nSD, nSpc
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

NSPC = 1
ISPC(1) = 1
iSM(1) = iSym
IAMCMP = 0
ICISTR = 1
i = 2
if (isym == state_sym) i = 1

if (NOCSF == 0) then
  nsd = max(ncsf(isym),nint(XISPSM(ISYM,1)))
  ipdcsfi = ipget(nsd)
  call ipin(ipdcsfi)
  ipDSDi = ipGet(nSD)
else
  nsd = max(ncsf(isym),nint(XISPSM(ISYM,1)))
  ipDSDi = ipGet(nsd)
  call ipin(ipdsdi)
end if

if (NOCSF == 0) then
  nD = NCSASM(ISYM)
  ipdiai = ipdcsfi
else
  nD = int(XISPSM(ISYM,ISPC(1)))
  ipdiai = ipdsdi
end if

call ipin(ipDSDi)
call IntDia(W(ipDSDi)%A,NSPC,ISPC,ISM,IAMCMP,rin_ene+potnuc)
if (NOCSF /= 1) call CSDIAG_MCLR(W(ipDCSFi)%A,W(ipDSDi)%A,NCNATS(1,ISYM),NTYP,CNSM(i)%ICTS,NDPCNT,NCPCNT)

if (NOCSF == 0) call ipclose(ipDSDi)

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
