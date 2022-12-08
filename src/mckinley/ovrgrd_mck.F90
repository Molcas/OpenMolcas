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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine OvrGrd_mck( &
#                     define _CALLING_
#                     include "grd_mck_interface.fh"
                     )
!***********************************************************************
!                                                                      *
! Object: to compute the gradients of the overlap matrix               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!     Author: Anders Bernhardsson                                      *
!             November '90                                             *
!                                                                      *
!             Modified to gradients of the overlap matrix. October     *
!             '91.                                                     *
!             Modified for respons calculation in May '95  By          *
!             Anders Bernhardsson                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp) :: iBeta, ip, ipAlph, ipAxyz, ipBeta, ipBxyz, ipRnxyz, ipRxyz, ipScrt, nip
logical(kind=iwp) :: ABeq(3)

#include "macros.fh"
unused_var(ZInv)
unused_var(lOper)
unused_var(iDCnt)
unused_var(iStabM)
unused_var(nStabM)

ABeq(1) = A(1) == RB(1)
ABeq(2) = A(2) == RB(2)
ABeq(3) = A(3) == RB(3)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipRnxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)*(nOrdOp+1)
ipAlph = nip
nip = nip+nZeta
ipBeta = nip
nip = nip+nZeta
ipScrt = nip
nip = nip+nTri_Elem1(la)*nTri_Elem1(lb)*nZeta*2

if (nip-1 > nArr) then
  write(u6,*) 'OvrGrd_Mck: nip-1 > nArr'
  write(u6,*) 'nip,nArr=',nip,nArr
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' IfGrad=',IfGrad
write(u6,*) ' IndGrd=',IndGrd
call RecPrt(' In OvrGrd_McK: A',' ',A,1,3)
call RecPrt(' In OvrGrd_McK: RB',' ',RB,1,3)
call RecPrt(' In OvrGrd_McK: Ccoor',' ',Ccoor,1,3)
call RecPrt(' In OvrGrd_McK: P',' ',P,nZeta,3)
write(u6,*) ' In OvrGrd_McK: la,lb=',la,lb
#endif

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipRnxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

! Combine the cartesian components to the gradient of the one
! electron integral and contract with the Fock matrix.

ip = ipAlph
do iBeta=1,nBeta
  Array(ip:ip+nAlpha-1) = Alpha
  ip = ip+nAlpha
end do
ip = ipBeta
do iBeta=1,nBeta
  Array(ip:ip+nAlpha-1) = Beta(iBeta)
  ip = ip+nAlpha
end do
call CmbnS1_mck(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,Array(ipScrt),Array(ipAlph),Array(ipBeta),IfGrad)

#ifdef _DEBUGPRINT_
call RecPrt(' Primitive Integrals',' ',Array(ipScrt),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
#endif

! Symmetry adapt the gradient operator

call SymAdO_mck(Array(ipScrt),nZeta*nTri_Elem1(la)*nTri_Elem1(lb),rFinal,nrOp,nop,IndGrd,iu,iv,ifgrad,idcar,trans)
#ifdef _DEBUGPRINT_
call RecPrt(' Primitive Integrals SO',' ',rFinal,nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nrOp)
#endif

return

end subroutine OvrGrd_mck
