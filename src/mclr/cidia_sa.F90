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

subroutine CIDIA_sa(iSym,ralp,S)

use Str_Info, only: CNSM
use ipPage, only: W
use MCLR_Data, only: ipCI
use MCLR_Data, only: ipDia
use MCLR_Data, only: FANCY_PRECONDITIONER
use MCLR_Data, only: XISPSM
use MCLR_Data, only: NOCSF, ICISTR
use MCLR_Data, only: NCNATS, NCPCNT, NDPCNT, NTYP
use input_mclr, only: State_Sym, rIn_Ene, PotNuc, ERASSCF, nCSF, nRoots, Weight
use Constants, only: Zero, One

implicit none
integer iSym
real*8 ralp(*), S(*)
integer iSM(1), LSPC(1), iSPC(1)
integer nSpc, iAMCmp, i, nSD, iPDCSFI, iPDSDI, ipDIAI, iP2, J
real*8 ECAS, WE
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

if (nocsf == 0) then
  ipdiai = ipdcsfi
else
  ipdiai = ipdsdi
end if
LSPC(1) = nSD

call ipin(ipDSDi)
call IntDia(W(ipDSDi)%Vec,NSPC,ISPC,ISM,LSPC,IAMCMP,rin_ene+potnuc)

if (NOCSF /= 1) call CSDIAG_MCLR(W(ipDCSFi)%Vec,W(ipDSDi)%Vec,NCNATS(1,ISYM),NTYP,CNSM(i)%ICTS,NDPCNT,NCPCNT)

if (nocsf == 0) call ipClose(ipDSDi)
! Calculate explicit part of hamiltonian

ipdia = ipdiai

if (FANCY_PRECONDITIONER) then
  call ipin(ipdia)
  call SA_PREC(S,W(ipdia)%Vec)
else
  call ipin(ipdiai)
  call ipin(ipCI)
  ip2 = 1
  do j=1,nroots
    ECAS = ERASSCF(j)
    We = Weight(j)
    ralp(j) = Zero
    do i=1,ncsf(State_SYM)
      ralp(j) = ralp(j)+One/(W(ipdiai)%Vec(i)-ECAS)*We*W(ipCI)%Vec(ip2)**2
      ip2 = ip2+1
    end do
  end do
end if

end subroutine CIDIA_sa
