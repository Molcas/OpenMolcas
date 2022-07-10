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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine MVeInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the mass-velocity integrals with the Gauss-       *
!         Hermite quadrature.                                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden. February '91                 *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, ib, iBeta, ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipQxyz, iPrint, iprV2, iprV4, ipRxyz, iRout, nip
logical(kind=iwp) :: ABeq(3)
character(len=80) :: Label

#include "macros.fh"
unused_var(ZInv)
unused_var(lOper)
unused_var(iChO)
unused_var(iStabM)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 190
iPrint = nPrint(iRout)
ABeq(:) = A == RB

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+3)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+3)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp-3)
ipQxyz = nip
nip = nip+nZeta*3*(la+3)*(lb+3)*(nOrdOp-3)
iprV2 = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*2
iprV4 = nip
nip = nip+nZeta*3*(la+1)*(lb+1)
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'MVeInt: nip-1 > nArr*nZeta')
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in MVeInt'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In MVeInt: A',' ',A,1,3)
  call RecPrt(' In MVeInt: RB',' ',RB,1,3)
  call RecPrt(' In MVeInt: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In MVeInt: P',' ',P,nZeta,3)
  call RecPrt(' In MVeInt: Zeta',' ',Zeta,nZeta,1)
  call RecPrt(' In MVeInt: Roots',' ',HerR(iHerR(nHer)),nHer,1)
  write(u6,*) ' In MVeInt: la,lb=',la,lb
end if

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+2,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+2,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(:) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp-4,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipQxyz),Array(ipAxyz),la+2,Array(ipRxyz),nOrdOp-4,Array(ipBxyz),lb+2,nZeta,HerW(iHerW(nHer)),nHer)

! Compute the cartesian components for the mass-velocity integrals.
! The kinetic energy components are linear combinations of overlap components.

ipAOff = ipA-1
do iBeta=1,nBeta
  Array(ipAOff+1:ipAOff+nAlpha) = Alpha
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB-1
do iBeta=1,nBeta
  Array(ipBOff+1:ipBOff+nAlpha) = Beta(iBeta)
  ipBOff = ipBOff+nAlpha
end do

call MVe(Array(iprV2),Array(iprV4),Array(ipQxyz),la,lb,Array(ipA),Array(ipB),nZeta)

! Combine the cartesian components to the full one electron integral.

call CmbnMV(Array(ipQxyz),nZeta,la,lb,nOrdOp-4,Zeta,rKappa,rFinal,nComp,Array(iprV2),Array(iprV4))

if (iPrint >= 99) then
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') 'Mass-Velocity(',ia,',',ib,')'
      call RecPrt(Label,' ',rFinal(:,ia,ib,:),nZeta,nComp)
    end do
  end do
end if

return

end subroutine MVeInt
