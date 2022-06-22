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

subroutine VeInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the velocity integrals with the Gauss-Hermite     *
!         quadrature.                                                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables.
logical ABeq(3)
integer iStabO(0:7), iDCRT(0:7)
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#include "macros.fh"
unused_var(Alpha)
unused_var(ZInv)
unused_var(nOrdOp)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 195
iPrint = nPrint(iRout)
ABeq(1) = A(1) == RB(1)
ABeq(2) = A(2) == RB(2)
ABeq(3) = A(3) == RB(3)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer
ipQxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+2)
ipVxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)
ipB = nip
nip = nip+nZeta
ipRes = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'VeInt: nip-1 > nArr*nZeta')
  write(6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(6,*) ' Abend in VeInt'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In VeInt: A',' ',A,1,3)
  call RecPrt(' In VeInt: RB',' ',RB,1,3)
  call RecPrt(' In VeInt: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In VeInt: P',' ',P,nZeta,3)
  write(6,*) ' In VeInt: la,lb=',la,lb
end if

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

! Compute the cartesian values of the basis functions angular part

call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the contribution from the multipole moment operator

ABeq(1) = .false.
ABeq(2) = .false.
ABeq(3) = .false.
call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),0,HerR(iHerR(nHer)),nHer,ABeq)

! Compute the cartesian components for the multipole moment
! integrals. The integrals are factorized into components.

call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),0,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

! Compute the cartesian components for the velocity integrals.
! The velocity components are linear combinations of overlap components.

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

call VelInt(Array(ipVxyz),Array(ipQxyz),la,lb,Array(ipB),nZeta)

! Combine the cartesian components to the full one electron integral.

call CmbnVe(Array(ipQxyz),nZeta,la,lb,0,Zeta,rKappa,Array(ipRes),nComp,Array(ipVxyz))

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine VeInt
