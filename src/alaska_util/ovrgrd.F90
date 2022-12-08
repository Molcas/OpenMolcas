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
!***********************************************************************

subroutine OvrGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the gradients of the overlap matrix               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to gradients of the overlap matrix. October     *
!             '91.                                                     *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Center_Info, only: dc
use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: iAlpha, iBeta, ip, ipAlph, ipAxyz, ipBeta, ipBxyz, iPrint, ipRnxyz, ipRxyz, iRout, nip
logical(kind=iwp) :: ABeq(3)
#include "print.fh"

#include "macros.fh"
unused_var(ZInv)
unused_var(iStabM)
unused_var(nStabM)

iRout = 122
iPrint = nPrint(iRout)
!write(u6,*) ' IfGrad=',IfGrad
!write(u6,*) ' IndGrd=',IndGrd
ABeq(:) = A == RB

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
if (nip-1 > nArr*nZeta) then
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in OvrGrd'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In OvrGrd: A',' ',A,1,3)
  call RecPrt(' In OvrGrd: RB',' ',RB,1,3)
  call RecPrt(' In OvrGrd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In OvrGrd: P',' ',P,nZeta,3)
  write(u6,*) ' In OvrGrd: la,lb=',la,lb
end if

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(:) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipRnxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

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
call CmbnS1(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,rFinal,Array(ipAlph),Array(ipBeta),Grad,nGrad,DAO,IfGrad,IndGrd,dc(mdc)%nStab, &
            dc(ndc)%nStab,kOp)

return

end subroutine OvrGrd
