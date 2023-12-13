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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine KnEInt_GIAO( &
#                      define _CALLING_
#                      include "int_interface.fh"
                      )
!***********************************************************************
!                                                                      *
! Object: to compute the kinetic energy integrals with the Gauss-      *
!         Hermite quadrature.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iBeta, iComp, iDCRT(0:7), ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipFnl, ipQxyz, iPrint, ipRxyz, ipTxyz, &
                     ipWxyz, iRout, iStabO(0:7), lDCRT, llOper, LmbdT, nB, nDCRT, nip, nOp, nStabO
real(kind=wp) :: TC(3)
logical(kind=iwp) :: ABeq(3)
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(ZInv)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 150
iPrint = nPrint(iRout)
ABeq(:) = A == RB

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+2)
ipQxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)*(nOrdOp+2)
ipTxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*(nOrdOp+2)
ipWxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*2
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
ipFnl = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
!                                                                      *
!***********************************************************************
!                                                                      *
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'KNEInt_GIAO: nip-1 > nArr*nZeta')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In KnEInt_GIAO: A',' ',A,1,3)
  call RecPrt(' In KnEInt_GIAO: RB',' ',RB,1,3)
  call RecPrt(' In KnEInt_GIAO: CoorO',' ',CoorO,1,3)
  call RecPrt(' In KnEInt_GIAO: P',' ',P,nZeta,3)
  write(u6,*) ' In KnEInt_GIAO: la,lb=',la,lb
end if
!                                                                      *
!***********************************************************************
!                                                                      *
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
!                                                                      *
!***********************************************************************
!                                                                      *

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CoorO,TC)

  ! Compute the contribution from the multipole moment operator

  ABeq(:) = .false.
  call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),nOrdOp+1,HerR(iHerR(nHer)),nHer,ABeq)

  ! Compute the cartesian components for the multipole moment
  ! integrals. The integrals are factorized into components.

  call Assmbl(Array(ipQxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp+1,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

  ! Compute the cartesian components for the kinetic energy
  ! integrals. The kinetic energy components are linear
  ! combinations of overlap components.

  ipAOff = ipA-1
  do iBeta=1,nBeta
    Array(ipAOff+1:ipAOff+nAlpha) = Alpha
    ipAOff = ipAOff+nAlpha
  end do

  ipBOff = ipB-1
  do iBeta=1,nBeta
    Array(ipBOff+1:ipBOff+nAlpha) = Beta(iBeta)
    ipBOff = ipBOff+nAlpha
  end do

  call Kntc_GIAO(Array(ipTxyz),Array(ipQxyz),Array(ipWxyz),la,lb,Array(ipA),Array(ipB),nZeta)

  ! Combine the cartesian components to the full one electron integral.

  nB = 3
  call CmbnKE_GIAO(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipFnl),nComp/nB,nB,Array(ipTxyz),Array(ipWxyz),A,RB,TC)

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine KnEInt_GIAO
