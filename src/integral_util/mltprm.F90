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
subroutine MltPrm( &
#                 define _CALLING_
#                 include "prm_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments integrals with the          *
!         Gauss-Hermite quadrature.                                    *
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
logical(kind=iwp) :: ABeq(3)
integer(kind=iwp) :: ipAxyz, ipBxyz, ipQxyz, ipRxyz, nip

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(iCnttp)

ABeq(:) = (A(:) == RB(:))

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipQxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'MltPrm: nip-1 > nArr*nZeta')
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in MltPrm'
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In MltPrm: A',' ',A,1,3)
call RecPrt(' In MltPrm: RB',' ',RB,1,3)
call RecPrt(' In MltPrm: Ccoor',' ',Ccoor,1,3)
call RecPrt(' In MltPrm: Kappa',' ',rKappa,nAlpha,nBeta)
call RecPrt(' In MltPrm: Zeta',' ',Zeta,nAlpha,nBeta)
call RecPrt(' In MltPrm: P',' ',P,nZeta,3)
write(u6,*) ' In MltPrm: la,lb=',la,lb
#endif

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(:) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

! Combine the cartesian components to the full one electron integral.

call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,rFinal,nComp)

return

end subroutine MltPrm
