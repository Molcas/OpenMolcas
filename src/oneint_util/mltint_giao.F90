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

subroutine MltInt_GIAO( &
#                      define _CALLING_
#                      include "int_interface.fh"
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

use Her_RW, only: HerR, HerW, iHerR, iHerW

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "oneswi.fh"
#include "print.fh"
#include "int_interface.fh"
real*8  RAB(3), TC(3)
character*80 Label, ChOper(0:7)*3
integer iStabO(0:7), iDCRT(0:7)
logical ABeq(3), EQ
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement function
nElem(i) = (i+1)*(i+2)/2

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 122
iPrint = nPrint(iRout)

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

if (.not. EQ(A,RB)) then

  ABeq(1) = A(1) == RB(1)
  ABeq(2) = A(2) == RB(2)
  ABeq(3) = A(3) == RB(3)
  RAB(1) = A(1)-RB(1)
  RAB(2) = A(2)-RB(2)
  RAB(3) = A(3)-RB(3)
  ! switch (only single center overlap matrix...)
  if (NDDO .and. (.not.(ABeq(1)) .and. ABeq(2) .and. ABeq(3))) then
    call dcopy_(nZeta*nIC*nElem(la)*nElem(lb),[Zero],0,final,1)
    return
  end if
  ! switch

  nip = 1
  ipAxyz = nip
  nip = nip+nZeta*3*nHer*(la+1)
  ipBxyz = nip
  nip = nip+nZeta*3*nHer*(lb+1)
  ipRxyz = nip
  nip = nip+nZeta*3*nHer*(nOrdOp+2)
  ipQxyz = nip
  nip = nip+nZeta*3*(la+1)*(lb+1)*(nOrdOp+2)
  ipFnl = nip
  nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
  if (nip-1 > nArr*nZeta) then
    call WarningMessage(2,'MltInt_GIAO: nip-1 > nArr*nZeta')
    write(6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
    write(6,*) ' Abend in MltInt'
    call Abend()
  end if

  if (iPrint >= 49) then
    call RecPrt(' In MltInt: A',' ',A,1,3)
    call RecPrt(' In MltInt: RB',' ',RB,1,3)
    call RecPrt(' In MltInt: Ccoor',' ',Ccoor,1,3)
    call RecPrt(' In MltInt: Kappa',' ',rKappa,nAlpha,nBeta)
    call RecPrt(' In MltInt: Zeta',' ',Zeta,nAlpha,nBeta)
    call RecPrt(' In MltInt: P',' ',P,nZeta,3)
    write(6,*) ' In MltInt: la,lb=',la,lb
  end if

  llOper = lOper(1)
  do iComp=2,nComp
    llOper = ior(llOper,lOper(iComp))
  end do

  ! Compute the cartesian values of the basis functions angular part

  call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
  call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb,HerR(iHerR(nHer)),nHer,ABeq)

  call SOS(iStabO,nStabO,llOper)
  call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
  if (iPrint >= 99) then
    write(6,*) ' m      =',nStabM
    write(6,'(9A)') '{M}=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(6,*) ' s      =',nStabO
    write(6,'(9A)') '{S}=',(ChOper(iStabO(ii)),ii=0,nStabO-1)
    write(6,*) ' LambdaT=',LmbdT
    write(6,*) ' t      =',nDCRT
    write(6,'(9A)') '{T}=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if

  do lDCRT=0,nDCRT-1
    call OA(iDCRT(lDCRT),CCoor,TC)

    ! Compute the contribution from the multipole moment operator

    ABeq(1) = .false.
    ABeq(2) = .false.
    ABeq(3) = .false.
    call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),nOrdOp+1,HerR(iHerR(nHer)),nHer,ABeq)

    ! Compute the cartesian components for the multipole moment
    ! integrals. The integrals are factorized into components.

    call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp+1,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

    ! Combine the cartesian components to the full one electron integral.

    nB = 3
    call CmbnMP_GIAO(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipFnl),nComp/nB,nB,RAB,TC)

    ! Accumulate contributions

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

  end do

end if

if (iPrint >= 99) then
  write(6,*) ' Result in MltInt'
  do ia=1,(la+1)*(la+2)/2
    do ib=1,(lb+1)*(lb+2)/2
      do iIC=1,nIC
        write(Label,'(A,I2,A,I2,A,I2,A)') ' Final(a=',ia,',b=',ib,',iIC=',iIC,')'
        call RecPrt(Label,' ',final(1,ia,ib,iIC),nAlpha,nBeta)
      end do
    end do
  end do
end if

return

end subroutine MltInt_GIAO
