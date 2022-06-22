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

subroutine MltInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
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
#include "rmat_option.fh"
#include "oneswi.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variable
real*8 TC(3), Origin(3)
integer iStabO(0:7), iDCRT(0:7)
logical ABeq(3), EQ
character*80 Label, ChOper(0:7)*3
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
data Origin/0.0d0,0.0d0,0.0d0/
! Statement function for Cartesian index
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

ABeq(1) = A(1) == RB(1)
ABeq(2) = A(2) == RB(2)
ABeq(3) = A(3) == RB(3)
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
nip = nip+nZeta*3*nHer*(nOrdOp+1)
ipQxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
ipFnl = nip
nip = nip+nZeta*nElem(la)*nElem(lb)*nComp
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_type_integrals) then
  ipRnr = nip
  nip = nip+nZeta*(la+lb+nOrdOp+1)
else
  ipRnr = -1
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'MltInt: nip-1 > nArr*nZeta')
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
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_type_integrals) then

  if (.not. EQ(CCoor,Origin)) then
    call WarningMessage(2,'MltInt: R-matrix error')
    write(6,*) 'MltInt: Wrong center of origin in case of R-matrix type of integrals!'
    write(6,*) ' Origin should always be (0.0,0.0,0.0)!'
    write(6,*) ' User the CENTER option to do this (see the SEWARD input section in the manual).'
    write(6,'(A,I3)') 'nOrdOp=',nOrdOp
    call Abend()
  end if

  ! R-matrix calculations: continuum basis functions (A=B=P=0)
  ! Compute the contributions of the basis functions and multipole
  ! radial part

  lsum = la+lb+nOrdOp
  call radlc(Zeta,nZeta,lsum,Array(ipRnr))

  ! Combine the radial and angular component to the full one electron integral.

  call CmbnMPr(Array(ipRnr),nZeta,la,lb,nOrdOp,Zeta,Array(ipFnl),nComp)

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

    ! Accumulate contributions

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)
  end do

else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
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
    call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

    ! Compute the cartesian components for the multipole moment
    ! integrals. The integrals are factorized into components.

    call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

    ! Combine the cartesian components to the full one electron integral.

    call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipFnl),nComp)

    ! Accumulate contributions

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

  end do

end if

if (iPrint >= 99) then
  write(6,*)
  write(6,*) ' Result in MltInt'
  write(6,*)
  write(6,*) 'la,lb,nHer=',la,lb,nHer
  write(6,*) 'nComp=',nComp
  write(6,*)
  do iIC=1,nIC
    write(Label,'(A,I2,A)') ' MltInt(iIC=',iIC,')'
    call RecPrt(Label,'(10G15.8) ',final(1,1,1,iIC),nZeta,nElem(la)*nElem(lb))
  end do
end if

return

end subroutine MltInt
