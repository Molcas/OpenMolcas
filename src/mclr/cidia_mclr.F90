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

use Exp, only: nexp, nexp_max
use Str_Info, only: CNSM
use ipPage, only: W
use negpre, only: nGP
use MCLR_Data, only: ipCI
use MCLR_Data, only: ipDia
use MCLR_Data, only: XISPSM
use MCLR_Data, only: NOCSF, ICISTR
use MCLR_Data, only: NCNATS, NCPCNT, NCSASM, NDPCNT, NTYP
use input_mclr, only: State_Sym, rIn_Ene, PotNuc, ERASSCF, nCSF, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp

implicit none
integer iSym
real*8 ralp
integer iSM(1), LSPC(1), iSPC(1), IDUM(1)
real*8, allocatable :: Q(:)
integer nSpc, iAMCmp, i, nSD, iPDCSFI, iRC, iPDSDI, nD, ipDIAI, nP2, nP1, nQ, iC
real*8 ECAS
real*8, external :: DDot_
integer, external :: ipClose, ipGet, ipIn, ipnOut

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
  irc = ipin(ipdcsfi)
  ipDSDi = ipGet(nSD)
else
  nsd = max(ncsf(isym),nint(XISPSM(ISYM,1)))
  ipDSDi = ipGet(nsd)
  irc = ipin(ipdsdi)
end if

if (NOCSF == 0) then
  nD = NCSASM(ISYM)
  ipdiai = ipdcsfi
else
  nD = idint(XISPSM(ISYM,ISPC(1)))
  ipdiai = ipdsdi
end if

LSPC(1) = nSD
irc = ipin(ipDSDi)
call IntDia(W(ipDSDi)%Vec,NSPC,ISPC,ISM,LSPC,IAMCMP,rin_ene+potnuc)
if (NOCSF /= 1) call CSDIAG_MCLR(W(ipDCSFi)%Vec,W(ipDSDi)%Vec,NCNATS(1,ISYM),NTYP,CNSM(i)%ICTS,NDPCNT,NCPCNT,0,0,IDUM)

if (NOCSF == 0) irc = ipclose(ipDSDi)

! Calculate explicit part of hamiltonian

np2 = min(nd,nexp_max)
np1 = 0
nq = 0
if (np2 /= 0) then
  irc = ipnout(ipdiai)
  irc = ipin(ipdiai)
  call h0(W(ipdiai)%Vec,np1,nexp_max,nq,isym,nexp,TimeDep)
else
  nexp = 0
end if

ECAS = ERASSCF(1)
irc = ipin(ipdiai)
do iC=1,nD
  if ((W(ipdiai)%Vec(ic)-ECAS) /= Zero) then
    W(ipdiai)%Vec(iC) = One/(W(ipdiai)%Vec(iC)-ECAS)
  else
    W(ipdiai)%Vec(iC) = 1.0e5_wp
  end if
end do
!             -1
! ralp=<0|(H-E) |0>

call mma_allocate(Q,nD,Label='Q')
Q(:) = Zero

irc = ipin(ipCI)
call ExpHinvv(W(ipdiai)%Vec,W(ipCI)%Vec,Q,Zero,One)

ralp = DDOT_(nD,W(ipCI)%Vec,1,Q,1)
if (NGP) then
  call MKP1INV(W(ipdiai)%Vec)
  call MKCIPRE()
end if
irc = ipnout(ipdiai)
call mma_deallocate(Q)

ipdia = ipdiai
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer(irc)
#endif

end subroutine CIDIA_MCLR
