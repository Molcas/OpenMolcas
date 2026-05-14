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

!#define _DEBUGPRINT_
subroutine KnEPrm( &
#                 define _CALLING_
#                 include "prm_interface.fh"
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

use Index_Functions, only: nTri_Elem1
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Definitions, only: wp, iwp, u6

implicit none
#include "prm_interface.fh"
integer(kind=iwp) :: nip, ipAxyz, ipBxyz, ipRxyz, ipQxyz, ipTxyz, ipA, ipB, ipAOff, ipBOff, iBeta
logical(kind=iwp) :: ABeq(3)

#include "macros.fh"
unused_var(ZInv)
unused_var(iCnttp)

ABeq(:) = (A(:) == RB(:))

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp-1)
ipQxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)*(nOrdOp-1)
ipTxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
!                                                                      *
!***********************************************************************
!                                                                      *
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'KnEPrm: nip-1 > nArr*nZeta')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In KnEPrm: A',' ',A,1,3)
call RecPrt(' In KnEPrm: RB',' ',RB,1,3)
call RecPrt(' In KnEPrm: Ccoor',' ',Ccoor,1,3)
call RecPrt(' In KnEPrm: P',' ',P,nZeta,3)
write(u6,*) ' In KnEPrm: la,lb=',la,lb
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp-2,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipQxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp-2,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

! Compute the cartesian components for the kinetic energy integrals.
! The kinetic energy components are linear combinations of overlap
! components.

ipAOff = ipA
ipBOff = ipB
do iBeta=1,nBeta
  Array(ipAOff:ipAOff+nAlpha-1) = Alpha(:)
  Array(ipBOff:ipBOff+nAlpha-1) = Beta(iBeta)
  ipAOff = ipAOff+nAlpha
  ipBOff = ipBOff+nAlpha
end do

call Kntc(Array(ipTxyz),Array(ipQxyz),la,lb,Array(ipA),Array(ipB),nZeta)

! Combine the cartesian components to the full one electron
! integral.

call CmbnKE(Array(ipQxyz),nZeta,la,lb,nOrdOp-2,Zeta,rKappa,rFinal,nComp,Array(ipTxyz))

return

end subroutine KnEPrm
