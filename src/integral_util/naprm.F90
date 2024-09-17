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

!#define _DEBUGPRINT_
subroutine NAPrm( &
#                define _CALLING_
#                include "prm_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January 1991                            *
!***********************************************************************

use Basis_Info, only: DBSC, Gaussian_Type, mGaussian_Type, Nuclear_Model, Point_Charge
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, OneHalf, Pi, TwoP54
use Definitions, only: wp, iwp

implicit none
#include "prm_interface.fh"
#include "oneswi.fh"
integer(kind=iwp) :: iAnga(4), ipIn, ipOff, lc, ld, mabMax, mabMin, mArr, mcdMax, mcdMin, nFlop, nMem, nT
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), EInv, Eta, Q_Nuc, rKappCD
logical(kind=iwp) :: NoSpecial
real(kind=wp), allocatable :: rKappa_mod(:)
logical(kind=iwp), external :: EQ
! Used for normal nuclear attraction integrals: TNAI, Fake, XCff2D, XRys2D
! Used for finite nuclei: TERI, ModU2, vCff2D, vRys2D
external :: Fake, ModU2, TERI, TNAI, vCff2D, vRys2D, XCff2D, XRys2D
! Statement function for Cartesian index
integer(kind=iwp) :: ixyz, nElem, nabSz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(nOrdOp)

call FZero(rFinal,nZeta*nElem(la)*nElem(lb)*nComp)

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
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1

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
  call mma_allocate(rKappa_mod,nZeta,label='rKappa_mod')
  rKappa_mod(:) = rKappa(:)*(TwoP54/Zeta(:))
end if

Q_Nuc = dbsc(iCnttp)%Charge

if (Q_Nuc /= Zero) then
  call dcopy_(3,CCoor,1,C,1)
# ifdef _DEBUGPRINT_
  call RecPrt('C',' ',C,1,3)
# endif

  call DCopy_(3,C,1,CoorAC(1,2),1)
  call DCopy_(3,C,1,Coori(1,3),1)
  call DCopy_(3,C,1,Coori(1,4),1)
  call DCopy_(3,C,1,Coora(1,3),1)
  call DCopy_(3,C,1,Coora(1,4),1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute integrals with the Rys quadrature.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nT = nZeta
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  select case (Nuclear_Model)
    case (Gaussian_Type)

      ! Gaussian nuclear charge distribution

      NoSpecial = .false.
      Eta = dbsc(iCnttp)%ExpNuc
      EInv = One/Eta
      rKappcd = TwoP54/Eta
      ! Tag on the normalization
      rKappcd = rKappcd*(Eta/Pi)**OneHalf
      ! s-type function
      mcdMin = 0
      mcdMax = 0
      call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa_mod,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin, &
               mcdMax,Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (mGaussian_Type)

      ! Modified Gaussian nuclear charge distribution

      NoSpecial = .false.
      Eta = dbsc(iCnttp)%ExpNuc
      EInv = One/Eta
      rKappcd = TwoP54/Eta
      ! Tag on the normalization
      rKappcd = rKappcd*(Eta/Pi)**OneHalf/(One+Three*dbsc(iCnttp)%w_mGauss/(Two*Eta))
      ! s type function
      mcdMin = 0
      mcdMax = 0
      call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa_mod,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin, &
               mcdMax,Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)

      ! d type function w*(x**2+y**2+z**2)
      if (dbsc(iCnttp)%w_mGauss > Zero) then
        rKappcd = rKappcd*dbsc(iCnttp)%w_mGauss
        iAnga(3) = 2
        mcdMin = nabSz(2+ld-1)+1
        mcdMax = nabSz(2+ld)
        ! tweak the pointers
        ipOff = 1+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
        mArr = nArr-(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
        call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin, &
                 mcdMax,Array(ipOff),mArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
        iAnga(3) = 0

        ! Add the s and d contributions together!

        call Assemble_mGauss(Array,Array(ipOff),nZeta*(mabMax-mabMin+1))
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    case (Point_Charge)

      ! Point-like nuclear charge distribution

      NoSpecial = .true.
      Eta = One
      EInv = One
      rKappcd = One
      mcdMin = 0
      mcdMax = 0
      call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin, &
               mcdMax,Array,nArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)
  end select
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Use the HRR to compute the required primitive integrals.

  call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

  call DCopy_(nZeta*nElem(la)*nElem(lb)*nComp,Array(ipIn),1,rFinal,1)
  call DScal_(nZeta*nElem(la)*nElem(lb)*nComp,-Q_Nuc,rFinal,1)
end if

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) call mma_deallocate(rKappa_mod)

return

end subroutine NAPrm
