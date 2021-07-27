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

subroutine KnEGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the kinetic energy integrals with the Gauss-      *
!         Hermite quadrature.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to gradients October '91.                       *
!***********************************************************************

use Her_RW
use Center_Info

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "grd_interface.fh"
! Local variables
logical ABeq(3)
! Statement function for Cartesian index
!Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
!iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
!Index(ixyz,ix,iz) = Ind(ixyz,ix,iz)+iOff(ixyz)

iRout = 150
iPrint = nPrint(iRout)
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
if (nip-1 > nArr*nZeta) then
  write(6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  call ErrTra()
  write(6,*) ' Abend in KnEGrd'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In KnEGrd: A',' ',A,1,3)
  call RecPrt(' In KnEGrd: RB',' ',RB,1,3)
  call RecPrt(' In KnEGrd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In KnEGrd: P',' ',P,nZeta,3)
  write(6,*) ' In KnEGrd: la,lb=',la,lb
end if

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
! The kinetic energy components are linear combinations of overlap
! components.

ipAOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

call Kntc(Array(ipTxyz),Array(ipRnxyz),la+1,lb+1,Array(ipA),Array(ipB),nZeta)

! Combine the cartesian components to the gradient of the kinetic
! energy integral and trace with the variational density matrix.

call CmbnT1(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,Final,Array(ipTxyz),Array(ipA),Array(ipB),Grad,nGrad,DAO,IfGrad,IndGrd, &
            dc(mdc)%nStab,dc(ndc)%nStab,kOp)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(ZInv)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iStabM)
end if

end subroutine KnEGrd
