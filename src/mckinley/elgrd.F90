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
! Copyright (C) 1990,1992,1995, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine ElGrd( &
#                define _CALLING_
#                include "1el_mck_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments integrals with the          *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!                                                                      *
!             Modified to reaction field calculations July '92         *
!             Modified to gradient calculations May '95                *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "1el_mck_interface.fh"
integer(kind=iwp) :: iBeta, ip, ipAlph, ipAxyz, ipBeta, ipBxyz, ipFinal, iprint, ipRnxyz, ipRxyz, ipTemp1, ipTemp2, ipTemp3, nip
logical(kind=iwp) :: ABeq(3)

iprint = 0
ABeq(:) = (A(:) == B(:))

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipRnxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)*(nOrdOp+1)
ipTemp1 = nip
nip = nip+nZeta
ipTemp2 = nip
nip = nip+nZeta
ipTemp3 = nip
nip = nip+3*nZeta*nHer
ipAlph = nip
nip = nip+nZeta
ipBeta = nip
nip = nip+nZeta
ipFinal = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*2
if (nip-1 > nArr) then
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr
  write(u6,*) ' Abend in ElGrd'
  call Abend()
end if

! Compute the cartesian values of the basis functions angular part

Array(ipTemp1:ipTemp2-1) = Zeta**(-Half)

call vCrtCmp(Array(ipTemp1),P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
call vCrtCmp(Array(ipTemp1),P,nZeta,B,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(:) = .false.
call vCrtCmp(Array(ipTemp1),P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call vAssmbl(Array(ipRnxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer,Array(ipTemp3))

! Combine the cartesian components to the full one electron integral.

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
call Cmbnel(Array(ipRnxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipFinal),Array(ipTemp1),Array(ipTemp2),Array(ipAlph), &
            Array(ipBeta),ifgrad,kcar)

!?
rFinal(:,:,:) = Zero

! Symmetry adapt the gradient operator

call SymAdO_mck2(Array(ipFinal),nZeta*nTri_Elem1(la)*nTri_Elem1(lb),rFinal,nrOp,nop,IndGrd,ksym,iu,iv,ifgrad,idcar,trans)
if (iPrint >= 49) call RecPrt(' Primitive Integrals SO',' ',rFinal,nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nrOp)

return

end subroutine ElGrd
