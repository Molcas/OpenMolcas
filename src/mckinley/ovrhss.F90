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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Ovrhss( &
#                 define _CALLING_
#                 include "hss_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the gradients of the overlap matrix               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Anders Bernhardsson 1995                                 *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Center_Info, only: dc
use Definitions, only: wp, iwp, u6

implicit none
#include "hss_interface.fh"
integer(kind=iwp) :: iAlpha, iBeta, ip, ipAlph, ipAxyz, ipBeta, ipBxyz, ipRnxyz, ipRxyz, nip
logical(kind=iwp) :: ABeq(3)

ABeq(1) = A(1) == RB(1)
ABeq(2) = A(2) == RB(2)
ABeq(3) = A(3) == RB(3)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+3)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+3)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipRnxyz = nip
nip = nip+nZeta*3*(la+3)*(lb+3)*(nOrdOp+1)
ipAlph = nip
nip = nip+nZeta
ipBeta = nip
nip = nip+nZeta
if (nip-1 > nArr) then
  write(u6,*) 'OvrHss: nip-1 > nArr'
  write(u6,*) 'nip,nArr=',nip,nArr
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' IfHss=',IfHss
write(u6,*) ' IndHss=',IndHss
call RecPrt(' In OvrHss: A',' ',A,1,3)
call RecPrt(' In OvrHss: RB',' ',RB,1,3)
call RecPrt(' In OvrHss: Ccoor',' ',Ccoor,1,3)
call RecPrt(' In OvrHss: P',' ',P,nZeta,3)
write(u6,*) ' In OvrHss: la,lb=',la,lb
#endif

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+2,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+2,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipRnxyz),Array(ipAxyz),la+2,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+2,nZeta,HerW(iHerW(nHer)),nHer)

! Combine the cartesian components to the gradient of the one
! electron integral and contract with the Fock matrix.

ip = ipAlph
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ip),1)
  ip = ip+nAlpha
end do
ip = ipBeta
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ip),nAlpha)
  ip = ip+1
end do

call CmbnS2(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,final,Array(ipAlph),Array(ipBeta),Hess,nHess,DAO,IfHss,IndHss, &
            indgrd,dc(mdc)%nStab,dc(ndc)%nStab,nOp)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(ZInv)
  call Unused_logical_array(ifgrd)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iStabM)
  call Unused_integer(nStabM)
end if

end subroutine Ovrhss
