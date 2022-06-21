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

subroutine NAInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January 1991                            *
!***********************************************************************

use Basis_Info
use Center_Info
use Gateway_global, only: Primitive_Pass
use DKH_Info, only: DKroll

implicit real*8(A-H,O-Z)
! Used for normal nuclear attraction integrals
external TNAI, Fake, XCff2D, XRys2D
! Used for finite nuclei
external TERI, ModU2, vCff2D, vRys2D
#include "real.fh"
#include "oneswi.fh"
#include "print.fh"
#include "int_interface.fh"
! Local arrys
real*8 C(3), TC(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
logical EQ, NoSpecial, No3Cnt, lECP
integer iAnga(4), iDCRT(0:7)
character ChOper(0:7)*3
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iRout = 151
iPrint = nPrint(iRout)

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

lECP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
end do
lc = 0
ld = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = lc
iAnga(4) = ld
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(2*3,Coora,1,Coori,1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
No3Cnt = .false.
if (EQ(A,RB)) then
  mabMin = nabSz(la+lb-1)+1
else if (NDDO) then
  No3Cnt = .true.
end if

! Compute FLOPs and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

! Modify Zeta if the two-electron code will be used!

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)*(TwoP54/Zeta(iZeta))
  end do
end if

! Loop over nuclear centers.

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp /= 1) kdc = kdc+dbsc(kCnttp-1)%nCntr

  ! Change nuclear charge if this is a relativistic ECP-case. This
  ! is used for the DKH transformation (see dkh_util/dkrelint_dp)!

  if (DKroll .and. Primitive_Pass .and. lECP) then
    Q_Nuc = dble(dbsc(kCnttp)%AtmNr)
  else
    Q_Nuc = dbsc(kCnttp)%Charge
  end if

  if (kCnttp == iCnttp_Dummy) cycle
  if (Q_Nuc == Zero) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
    if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    if (iPrint >= 99) then
      write(6,*) ' m      =',nStabM
      write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
      write(6,*) ' s      =',dc(kdc+kCnt)%nStab
      write(6,'(9A)') '(S)=',(ChOper(dc(kdc+kCnt)%iStab(ii)),ii=0,dc(kdc+kCnt)%nStab-1)
      write(6,*) ' LambdaT=',LmbdT
      write(6,*) ' t      =',nDCRT
      write(6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
    end if

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      ! switch (only two center NA matrix...)
      if (No3Cnt .and. (.not. (EQ(A,TC) .or. EQ(RB,TC)))) Go To 102
      ! switch
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      call dcopy_(3,TC,1,Coora(1,3),1)
      call dcopy_(3,TC,1,Coora(1,4),1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      !        Compute integrals with the Rys quadrature.
      !                                                                *
      !*****************************************************************
      !                                                                *
      nT = nZeta
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (Nuclear_Model == Gaussian_Type) then

        ! Gaussian nuclear charge distribution

        NoSpecial = .false.
        Eta = dbsc(kCnttp)%ExpNuc
        EInv = One/Eta
        rKappcd = TwoP54/Eta
        ! Tag on the normalization
        rKappcd = rKappcd*(Eta/Pi)**(Three/Two)
        ! s-type function
        mcdMin = 0
        mcdMax = 0
        call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin, &
                 mcdMax,Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
        !                                                              *
        !***************************************************************
        !                                                              *
      else if (Nuclear_Model == mGaussian_Type) then

        ! Modified Gaussian nuclear charge distribution

        NoSpecial = .false.
        Eta = dbsc(kCnttp)%ExpNuc
        EInv = One/Eta
        rKappcd = TwoP54/Eta
        ! Tag on the normalization
        rKappcd = rKappcd*(Eta/Pi)**(Three/Two)/(One+Three*dbsc(kCnttp)%w_mGauss/(Two*Eta))
        ! s type function
        mcdMin = 0
        mcdMax = 0
        call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin, &
                 mcdMax,Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)

        ! d type function w*(x**2+y**2+z**2)
        if (dbsc(kCnttp)%w_mGauss > 0.0d0) then
          rKappcd = rKappcd*dbsc(kCnttp)%w_mGauss
          iAnga(3) = 2
          mcdMin = nabSz(2+ld-1)+1
          mcdMax = nabSz(2+ld)
          ! tweak the pointers
          ipOff = 1+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
          mArr = nArr-(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
          call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin, &
                   mcdMax,Array(ipOff),mArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
          iAnga(3) = 0

          ! Add the s and d contributions together!

          call Assemble_mGauss(Array,Array(ipOff),nZeta*(mabMax-mabMin+1))
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
      else if (Nuclear_Model == Point_Charge) then

        ! Point-like nuclear charge distribution

        NoSpecial = .true.
        Eta = One
        EInv = One
        rKappcd = One
        mcdMin = 0
        mcdMax = 0
        call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin, &
                 mcdMax,Array,nArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      !
      ! Use the HRR to compute the required primitive integrals.

      call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

      ! Accumulate contributions to the symmetry adapted operator

      nOp = NrOpr(iDCRT(lDCRT))
      call SymAdO(Array(ipIn),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,-Fact*Q_Nuc)
      if (iPrint >= 99) then
        write(6,*) Fact*Q_Nuc
        call RecPrt('NaInt: Array(ipIn)',' ',Array(ipIn),nZeta,nElem(la)*nElem(lb)*nComp)
        call RecPrt('NaInt: Final',' ',final,nZeta,nElem(la)*nElem(lb)*nIC)
      end if

102   continue
    end do
  end do
end do

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)/(TwoP54/Zeta(iZeta))
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_integer(nHer)
  call Unused_real_array(CCoor)
  call Unused_integer(nOrdOp)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine NAInt
