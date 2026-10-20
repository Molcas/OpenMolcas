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
! Copyright (C) 1993, Bernd Artur Hess                                 *
!***********************************************************************

subroutine VPInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of  pV integrals          *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, April 1993          *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Integral_interfaces, only: int_kernel
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: iBeta, ipArr, ipB, ipOff, ipS1, ipS2, kComp, kIC, kRys, mArr, nip, nRys
procedure(int_kernel) :: NAint
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i

call RecPrt(' In VpInt: Alpha','(5ES20.13)',Alpha,nAlpha,1)
call RecPrt(' In VpInt: Beta','(5ES20.13)',Beta,nBeta,1)
#endif

nRys = nHer

nip = 1
ipB = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb+1)
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb-1)
else
  ipS2 = ipS1
end if
ipArr = nip
mArr = nArr-(nip-1)/nZeta
if (mArr < 0) then
  call WarningMessage(2,'VpInt: mArr<0!')
  call Abend()
end if

rFinal(:,:,:,:) = Zero
Array(:) = Zero
! Compute contribution from a,b+1

kRys = ((la+1)+lb+2)/2

kIC = 1
kComp = 1
call NAInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,nIC,nComp,la,lb+1,A,RB,kRys,Array(ipArr),mArr,CoorO, &
           nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)

ipOff = ipB-1
do iBeta=1,nBeta
  Array(ipOff+1:ipOff+nAlpha) = Beta(iBeta)
  ipOff = ipOff+nAlpha
end do

! Compute contribution from a,b-1

if (lb > 0) then
  kRys = ((la-1)+lb+2)/2

  call NAInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,nRys,Array(ipArr),mArr,CoorO, &
             nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
end if

! Assemble final integral from the derivative integrals

#ifdef _DEBUGPRINT_
call RecPrt(' In VpInt: Beta (expanded)','(5ES20.13)',Array(ipB),nZeta,1)
#endif

call Util8(Array(ipB),nZeta,rFinal,la,lb,Array(ipS1),Array(ipS2))

#ifdef _DEBUGPRINT_
do i=1,3
  call RecPrt('VpInt: rFinal',' ',rFinal(:,:,:,i),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
end do
#endif

end subroutine VPInt
