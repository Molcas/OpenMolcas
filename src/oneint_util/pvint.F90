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

subroutine PVInt( &
#                define _CALLING_
#                include "int_interface.fh"
                 ,Kernel)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of  pX integrals          *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, April 1993          *
!***********************************************************************

implicit real*8(A-H,O-Z)
external Kernel
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Statement function for Cartesian index
nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 221
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*) 'PVInt: nIC,nComp=',nIC,nComp
  call RecPrt(' In pvint: Alpha','(5D20.13)',Alpha,nAlpha,1)
  call RecPrt(' In pvint: Beta','(5D20.13)',Beta,nBeta,1)
end if

nip = 1
ipA = nip
nip = nip+nZeta
ipS1 = nip
nip = nip+nZeta*nElem(la+1)*nElem(lb)*nIC
ipS2 = 1
if (la > 0) then
  ipS2 = nip
  nip = nip+nZeta*nElem(la-1)*nElem(lb)*nIC
else
  ipS2 = ipS1
end if
ipArr = nip
mArr = nArr-(nip-1)/nZeta
if (mArr < 0) then
  call WarningMessage(2,'pVInt: mArr<0!')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute contribution from a+1,b

kRys = ((la+1)+lb+2)/2
call Kernel(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS1),nZeta,nIC,nComp,la+1,lb,A,RB,kRys,Array(ipArr),mArr,CCoor, &
            nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute contribution from a-1,b

if (la > 0) then
  kRys = ((la-1)+lb+2)/2
  call Kernel(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,Array(ipS2),nZeta,nIC,nComp,la-1,lb,A,RB,kRys,Array(ipArr),mArr,CCoor, &
              nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
ipOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipOff),1)
  ipOff = ipOff+nAlpha
end do
if (iPrint >= 99) then
  call RecPrt(' In pvint: Alpha (expanded)','(5D20.13)',Array(ipA),nZeta,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Assemble final integral from the derivative integrals

call Ass_pX(Array(ipA),nZeta,final,la,lb,Array(ipS1),Array(ipS2),nIC)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 49) then
  do i=1,3
    call RecPrt('pVInt: Final',' ',final(1,1,1,i),nZeta,nElem(la)*nElem(lb))
  end do
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nHer)

end subroutine PVInt
