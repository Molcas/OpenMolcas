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

use Basis_Info, only: dbsc, Gaussian_Type, iCnttp_Dummy, mGaussian_Type, nCnttp, Nuclear_Model, Point_Charge
use Center_Info, only: dc
use Gateway_global, only: Primitive_Pass
use DKH_Info, only: DKroll
use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Constants, only: Zero, One, Two, Three, OneHalf, Pi, TwoP54
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "oneswi.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAnga(4), iDCRT(0:7), ii, ipIn, ipOff, iPrint, iRout, kCnt, kCnttp, kdc, lc, ld, lDCRT, LmbdT, mabMax, &
                     mabMin, mArr, mcdMax, mcdMin, nDCRT, nFLOP, nMem, nOp, nT
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), EInv, Eta, Fact, Q_Nuc, rKappcd, TC(3)
logical(kind=iwp) :: lECP, No3Cnt, NoSpecial
character(len=*), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
! Used for normal nuclear attraction integrals: TNAI, Fake, XCff2D, XRys2D
! Used for finite nuclei: TERI, ModU2, vCff2D, vRys2D
external :: Fake, ModU2, TERI, TNAI, vCff2D, vRys2D, XCff2D, XRys2D

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 151
iPrint = nPrint(iRout)

rFinal(:,:,:,:) = Zero

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
Coora(:,1) = A
Coora(:,2) = RB
Coori(:,1:2) = Coora(:,1:2)
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
No3Cnt = .false.
if (EQ(A,RB)) then
  mabMin = nTri3_Elem1(la+lb-1)
else if (NDDO) then
  No3Cnt = .true.
end if

! Compute FLOPs and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if

! Modify Zeta if the two-electron code will be used!

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  rKappa(:) = rKappa*(TwoP54/Zeta)
end if

! Loop over nuclear centers.

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr

  ! Change nuclear charge if this is a relativistic ECP-case. This
  ! is used for the DKH transformation (see dkh_util/dkrelint_dp)!

  if (DKroll .and. Primitive_Pass .and. lECP) then
    Q_Nuc = real(dbsc(kCnttp)%AtmNr,kind=wp)
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
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    if (iPrint >= 99) then
      write(u6,*) ' m      =',nStabM
      write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
      write(u6,*) ' s      =',dc(kdc+kCnt)%nStab
      write(u6,'(9A)') '(S)=',(ChOper(dc(kdc+kCnt)%iStab(ii)),ii=0,dc(kdc+kCnt)%nStab-1)
      write(u6,*) ' LambdaT=',LmbdT
      write(u6,*) ' t      =',nDCRT
      write(u6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
    end if

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      ! switch (only two center NA matrix...)
      if (No3Cnt .and. (.not. (EQ(A,TC) .or. EQ(RB,TC)))) cycle
      ! switch
      CoorAC(:,2) = TC
      Coori(:,3) = TC
      Coori(:,4) = TC
      Coora(:,3) = TC
      Coora(:,4) = TC
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
        rKappcd = rKappcd*(Eta/Pi)**OneHalf
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
        rKappcd = rKappcd*(Eta/Pi)**OneHalf/(One+Three*dbsc(kCnttp)%w_mGauss/(Two*Eta))
        ! s type function
        mcdMin = 0
        mcdMax = 0
        call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin, &
                 mcdMax,Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)

        ! d type function w*(x**2+y**2+z**2)
        if (dbsc(kCnttp)%w_mGauss > Zero) then
          rKappcd = rKappcd*dbsc(kCnttp)%w_mGauss
          iAnga(3) = 2
          mcdMin = nTri3_Elem1(2+ld-1)
          mcdMax = nTri3_Elem1(2+ld)-1
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
      call SymAdO(Array(ipIn),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,-Fact*Q_Nuc)
      if (iPrint >= 99) then
        write(u6,*) Fact*Q_Nuc
        call RecPrt('NaInt: Array(ipIn)',' ',Array(ipIn),nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nComp)
        call RecPrt('NaInt: rFinal',' ',rFinal,nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nIC)
      end if

    end do
  end do
end do

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  rKappa = rKappa/(TwoP54/Zeta)
end if

return

end subroutine NAInt
