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
use ipPage, only: ipclose, ipget, ipin, W
use MCLR_Data, only: FANCY_PRECONDITIONER, ICISTR, ipCI, ipDia, NCNATS, NCPCNT, NDPCNT, NOCSF, NTYP, XISPSM
use input_mclr, only: ERASSCF, nCSF, nRoots, PotNuc, rIn_Ene, State_Sym, Weight
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(inout) :: ralp(nRoots), S(*)
integer(kind=iwp) :: i, iAMCmp, iP2, iPDCSFI, ipDIAI, iPDSDI, iSM(1), iSPC(1), J, nSD, nSpc

! This is just a interface to hide Jeppe from the rest of the world
! we don't want to let world see the work of the Danish
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

call ipin(ipDSDi)
call IntDia(W(ipDSDi)%A,NSPC,ISPC,ISM,IAMCMP,rin_ene+potnuc)

if (NOCSF /= 1) call CSDIAG_MCLR(W(ipDCSFi)%A,W(ipDSDi)%A,NCNATS(1,ISYM),NTYP,CNSM(i)%ICTS,NDPCNT,NCPCNT)

if (nocsf == 0) call ipClose(ipDSDi)
! Calculate explicit part of hamiltonian

ipdia = ipdiai

if (FANCY_PRECONDITIONER) then
  call ipin(ipdia)
  call SA_PREC(S,W(ipdia)%A)
else
  call ipin(ipdiai)
  call ipin(ipCI)
  ip2 = 0
  do j=1,nroots
    ralp(j) = sum(One/(W(ipdiai)%A(1:ncsf(State_Sym))-ERASSCF(j))*Weight(j)*W(ipCI)%A(ip2+1:ip2+ncsf(State_Sym))**2)
    ip2 = ip2+ncsf(State_Sym)
  end do
end if

end subroutine CIDIA_sa
