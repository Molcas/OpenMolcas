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

implicit real*8(A-H,O-Z)
external TNAI, Fake, XCff2D, XRys2D
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Statement function for Cartesian index
nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2

iRout = 221
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In vpint: Alpha','(5D20.13)',Alpha,nAlpha,1)
  call RecPrt(' In vpint: Beta','(5D20.13)',Beta,nBeta,1)
end if

nRys = nHer

nip = 1
ipB = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nElem(la)*nElem(lb+1)
if (lb > 0) then
  ipS2 = nip
  nip = nip+nZeta*nElem(la)*nElem(lb-1)
else
  ipS2 = ipS1
end if
ipArr = nip
mArr = nArr-(nip-1)/nZeta
if (mArr < 0) then
  call WarningMessage(2,'VpInt: mArr<0!')
  call Abend()
end if

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)
call dcopy_(nZeta*nArr,[Zero],0,Array,1)
! Compute contribution from a,b+1

kRys = ((la+1)+lb+2)/2

kIC = 1
kComp = 1
call NAInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,nIC,nComp,la,lb+1,A,RB,kRys,Array(ipArr),mArr,CCoor, &
           nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)

ipOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
  ipOff = ipOff+1
end do

! Compute contribution from a,b-1

if (lb > 0) then
  kRys = ((la-1)+lb+2)/2

  call NAInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,nRys,Array(ipArr),mArr,CCoor, &
             nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
end if

! Assemble final integral from the derivative integrals

if (iPrint >= 99) call RecPrt(' In vpint: Beta (expanded)','(5D20.13)',Array(ipB),nZeta,1)

call Util8(Array(ipB),nZeta,final,la,lb,Array(ipS1),Array(ipS2))

if (iPrint >= 49) then
  do i=1,3
    call RecPrt('VpInt: Final',' ',final(1,1,1,i),nZeta,nElem(la)*nElem(lb))
  end do
end if

return

end subroutine VPInt
