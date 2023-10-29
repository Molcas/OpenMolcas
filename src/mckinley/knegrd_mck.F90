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

subroutine KnEGrd_mck( &
#                     define _CALLING_
#                     include "grd_mck_interface.fh"
                     )
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the kinetic energy integrals      *
!         with the Gauss-Hermite quadrature                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Anders Bernhardsson,1995                                 *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp) :: iBeta, ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipRnxyz, ipRxyz, ipSc, ipTxyz, nip
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
nip = nip+nZeta*3*nHer*(la+3)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+3)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipRnxyz = nip
nip = nip+nZeta*3*(la+3)*(lb+3)*(nOrdOp+1)
ipTxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
ipSc = nip
nip = nip+nTri_Elem1(la)*nTri_Elem1(lb)*nZeta
if (nip-1 > nArr) then
  write(u6,*) 'KneGrd_Mck: nip-1 > nArr'
  write(u6,*) 'nip,nArr=',nip,nArr
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In KnEGrd_McK: A',' ',A,1,3)
call RecPrt(' In KnEGrd_McK: RB',' ',RB,1,3)
call RecPrt(' In KnEGrd_McK: Ccoor',' ',Ccoor,1,3)
call RecPrt(' In KnEGrd_McK: P',' ',P,nZeta,3)
write(u6,*) ' In KnEGrd_McK: la,lb=',la,lb
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

! Compute the cartesian components for the kinetic energy integrals.
! The kinetic energy components are linear combinations of overlap components.

ipAOff = ipA
do iBeta=1,nBeta
  Array(ipAOff:ipAOff+nAlpha-1) = Alpha
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iBeta=1,nBeta
  Array(ipBOff:ipBOff+nAlpha-1) = Beta(iBeta)
  ipBOff = ipBOff+nAlpha
end do

call Kntc(Array(ipTxyz),Array(ipRnxyz),la+1,lb+1,Array(ipA),Array(ipB),nZeta)

! Combine the cartesian components to the gradient of the kinetic
! energy integral and trace with the variational density matrix.

call CmbnT1_mck(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,Array(ipSc),Array(ipTxyz),Array(ipA),Array(ipB),IfGrad)

rFinal(:,:,:,:) = Zero

! Symmetry adapt the gradient operator

call SymAdO_mck(Array(ipSc),nZeta*nTri_Elem1(la)*nTri_Elem1(lb),rFinal,nrOp,nop,IndGrd,iu,iv,ifgrad,idCar,trans)

return

end subroutine KnEGrd_mck
