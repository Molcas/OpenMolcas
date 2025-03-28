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

subroutine CIDIA_TD(iSym)

use Exp, only: nexp, nexp_max
use Str_Info, only: CNSM
use ipPage, only: W
use MCLR_Data, only: ipDia
use MCLR_Data, only: XISPSM
use MCLR_Data, only: NOCSF, ICISTR
use MCLR_Data, only: NCNATS, NCPCNT, NCSASM, NDPCNT, NTYP
use input_mclr, only: State_Sym, rIn_Ene, PotNuc, ERASSCF, nCSF, TimeDep

implicit none
integer iSym
integer iSM(1), LSPC(1), iSPC(1)
integer nSpc, iAMCmp, i, nSD, iPDCSFI, iPDSDI, nD, ipDIAI, nP2, nP1, nQ, iC
real*8 ECAS
integer, external :: ipGet

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
  nD = idint(XISPSM(ISYM,ISPC(1)))
  ipdiai = ipdsdi
end if

LSPC(1) = nSD

call ipin(ipDSDi)
call IntDia(W(ipDSDi)%Vec,NSPC,ISPC,ISM,LSPC,IAMCMP,rin_ene+potnuc)
if (NOCSF /= 1) call CSDIAG_MCLR(W(ipDCSFi)%Vec,W(ipDSDi)%Vec,NCNATS(1,ISYM),NTYP,CNSM(i)%ICTS,NDPCNT,NCPCNT)

if (nocsf == 0) call ipClose(ipDSDi)
! Calculate explicit part of hamiltonian

np2 = min(nd,nexp_max)
np1 = 0
nq = 0
if (np2 /= 0) then
  call ipnout(ipdiai)
  call ipin(ipdiai)
  call h0(W(ipdiai)%Vec,np1,nexp_max,nq,isym,nexp,TimeDep)
else
  nexp = 0
end if

ECAS = ERASSCF(1)
call ipin(ipdiai)
do iC=1,nD
  W(ipdiai)%Vec(iC) = (W(ipdiai)%Vec(iC)-ECAS)
end do

ipdia = ipdiai

end subroutine CIDIA_TD
