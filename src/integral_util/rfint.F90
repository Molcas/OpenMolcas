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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine RFInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,nComp,la,lb,A,B,nHer,Array,nArr,Ccoor,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments integrals with the          *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to reaction field calculations July '92         *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerw
use Definitions, only: u6

implicit none
integer nAlpha, nBeta, nZeta, nComp, la, lb, nHer, nArr, nOrdOp
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), B(3), Array(nZeta*nArr), Ccoor(3)
logical ABeq(3)
integer nip, ipAxyz, ipBxyz, ipRnxyz, ipTemp1, ipTemp2, ipTemp3, iZeta, ipRxyz

ABeq(:) = A(:) == B(:)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipRnxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
ipTemp1 = nip
nip = nip+nZeta
ipTemp2 = nip
nip = nip+nZeta
ipTemp3 = nip
nip = nip+3*nZeta*nHer
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'RFInt: nip-1 > nArr*nZeta')
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in RFInt'
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In RFInt: A',' ',A,1,3)
call RecPrt(' In RFInt: B',' ',B,1,3)
call RecPrt(' In RFInt: CCoor',' ',CCoor,1,3)
call RecPrt(' In RFInt: P',' ',P,nZeta,3)
write(u6,*) ' In RFInt: la,lb=',la,lb
write(u6,*) ' In RFInt: nHer=',nHer
#endif

! Compute the cartesian values of the basis functions angular part

do iZeta=1,nZeta
  Array(ipTemp1-1+iZeta) = 1/sqrt(Zeta(iZeta))
  !Array(ipTemp1-1+iZeta) = Zeta(iZeta)**(-Half)
end do
call vCrtCmp(Array(ipTemp1),P,nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
call vCrtCmp(Array(ipTemp1),P,nZeta,B,Array(ipBxyz),lb,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(:) = .false.
call vCrtCmp(Array(ipTemp1),P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call vAssmbl(Array(ipRnxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer,Array(ipTemp3))

! Combine the cartesian components to the full one electron integral.

call CmbnRF(Array(ipRnxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,final,nComp,Array(ipTemp1),Array(ipTemp2))

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_real_array(ZInv)
end if
#endif

end subroutine RFInt
