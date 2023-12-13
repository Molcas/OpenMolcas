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
use Index_Functions, only: nTri_Elem1
use rmat, only: RMat_Type_Integrals
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "oneswi.fh"
#include "print.fh"
integer(kind=iwp) :: iComp, iDCRT(0:7), ii, iIC, ipAxyz, ipBxyz, ipFnl, ipQxyz, iPrint, ipRnr, ipRxyz, iRout, iStabO(0:7), lDCRT, &
                     llOper, LmbdT, lsum, nDCRT, nip, nOp, nStabO
real(kind=wp) :: TC(3)
logical(kind=iwp) :: ABeq(3)
character(len=80) :: Label
real(kind=wp), parameter :: Origin(3) = Zero
character(len=*), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 122
iPrint = nPrint(iRout)

rFinal(:,:,:,:) = Zero

ABeq(:) = A == RB
! switch (only single center overlap matrix...)
if (NDDO .and. (.not.(ABeq(1)) .and. ABeq(2) .and. ABeq(3))) return
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
nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
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
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in MltInt'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In MltInt: A',' ',A,1,3)
  call RecPrt(' In MltInt: RB',' ',RB,1,3)
  call RecPrt(' In MltInt: CoorO',' ',CoorO,1,3)
  call RecPrt(' In MltInt: Kappa',' ',rKappa,nAlpha,nBeta)
  call RecPrt(' In MltInt: Zeta',' ',Zeta,nAlpha,nBeta)
  call RecPrt(' In MltInt: P',' ',P,nZeta,3)
  write(u6,*) ' In MltInt: la,lb=',la,lb
end if

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_type_integrals) then

  if (.not. EQ(CoorO,Origin)) then
    call WarningMessage(2,'MltInt: R-matrix error')
    write(u6,*) 'MltInt: Wrong center of origin in case of R-matrix type of integrals!'
    write(u6,*) ' Origin should always be (0.0,0.0,0.0)!'
    write(u6,*) ' User the CENTER option to do this (see the SEWARD input section in the manual).'
    write(u6,'(A,I3)') 'nOrdOp=',nOrdOp
    call Abend()
  end if

  ! R-matrix calculations: continuum basis functions (A=B=P=0)
  ! Compute the contributions of the basis functions and multipole
  ! radial part

  lsum = la+lb+nOrdOp
  call radlc(Zeta,nZeta,lsum,Array(ipRnr))

  ! Combine the radial and angular component to the full one electron integral.

  call CmbnMPr(Array(ipRnr),nZeta,la,lb,nOrdOp,Array(ipFnl),nComp)

  call SOS(iStabO,nStabO,llOper)
  call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
  if (iPrint >= 99) then
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '{M}=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(u6,*) ' s      =',nStabO
    write(u6,'(9A)') '{S}=',(ChOper(iStabO(ii)),ii=0,nStabO-1)
    write(u6,*) ' LambdaT=',LmbdT
    write(u6,*) ' t      =',nDCRT
    write(u6,'(9A)') '{T}=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if

  do lDCRT=0,nDCRT-1

    ! Accumulate contributions

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)
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
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '{M}=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(u6,*) ' s      =',nStabO
    write(u6,'(9A)') '{S}=',(ChOper(iStabO(ii)),ii=0,nStabO-1)
    write(u6,*) ' LambdaT=',LmbdT
    write(u6,*) ' t      =',nDCRT
    write(u6,'(9A)') '{T}=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if

  do lDCRT=0,nDCRT-1
    call OA(iDCRT(lDCRT),CoorO,TC)

    ! Compute the contribution from the multipole moment operator

    ABeq(:) = .false.
    call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

    ! Compute the cartesian components for the multipole moment
    ! integrals. The integrals are factorized into components.

    call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

    ! Combine the cartesian components to the full one electron integral.

    call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipFnl),nComp)

    ! Accumulate contributions

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

  end do

end if

if (iPrint >= 99) then
  write(u6,*)
  write(u6,*) ' Result in MltInt'
  write(u6,*)
  write(u6,*) 'la,lb,nHer=',la,lb,nHer
  write(u6,*) 'nComp=',nComp
  write(u6,*)
  do iIC=1,nIC
    write(Label,'(A,I2,A)') ' MltInt(iIC=',iIC,')'
    call RecPrt(Label,'(10G15.8) ',rFinal(:,:,:,iIC),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
  end do
end if

return

end subroutine MltInt
