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
! Copyright (C) 1991,1995, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
subroutine EFInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electric field         *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
! Modified for explicit code, R. Lindh, February '95.                  *
!***********************************************************************

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: i, iAnga(4), iComp, iDCRT(0:7), ip1, ip2, ip3, ipIn, iStabO(0:7), kab, lab, labcd, lcd, lDCRT, llOper, LmbdT, &
                     mabMax, mabMin, mArr, mcdMax, mcdMin, nDCRT, nFLOP, nMem, nOp, nStabO, nT, nzab
real(kind=wp) :: CoorAC(3,2), Coori(3,4), RR, TC(3), XX, YY
logical(kind=iwp) :: NoSpecial
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iElem, ij, ip, jElem
character(len=80) :: Label
#endif
real(kind=wp), pointer :: EFInts(:,:)
real(kind=wp), parameter :: ThreeI = One/Three
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
external :: Fake, TNAI, XCff2D, XRys2D

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(PtChrg)
unused_var(iAddPot)

rFinal(:,:,:,:) = Zero

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
Coori(:,1) = A
Coori(:,2) = RB
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
mcdMin = nTri3_Elem1(nOrdOp-1)
mcdMax = nTri3_Elem1(nOrdop)-1
lab = (mabMax-mabMin+1)
kab = nTri_Elem1(la)*nTri_Elem1(lb)
lcd = (mcdMax-mcdMin+1)
labcd = lab*lcd

! Compute Flop's and size of work array which HRR will Use.

call mHRR(la,lb,nFLOP,nMem)

! Distribute the work array

ip2 = 1
ip1 = ip2+nZeta*max(labcd,lcd*nMem)
mArr = nArr-max(labcd,lcd*nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CoorO,TC)
  CoorAC(:,2) = TC(:)
  Coori(:,3) = TC(:)
  Coori(:,4) = TC(:)

  ! Compute integrals with the Rys-Gauss quadrature.

  nT = nZeta
  NoSpecial = .true.
  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
           Array(ip1),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

  ! The integrals are now ordered as ijkl,e,f

  ! a) Change the order to f,ijkl,e
  ! b) Unfold e to ab, f,ijkl,ab
  ! c) Change the order back to ijkl,ab,f

  ! a)

  call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)

  ! b) Use the HRR to unfold e to ab

  call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
  ip3 = ip2-1+ipIn

  ! c)

  call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,Array(ip1),nZeta*kab)

  ! Modify to traceless form, the sixth element contains r*r and

  if (nOrdOp == 2) then
    nzab = nZeta*kab
    EFInts(1:nzab,1:6) => Array(ip1:ip1-1+6*nzab)
    do i=1,nzab
      RR = EFInts(i,1)+EFInts(i,4)+EFInts(i,6)
      XX = Two*EFInts(i,1)-EFInts(i,4)-EFInts(i,6)
      YY = -EFInts(i,1)+Two*EFInts(i,4)-EFInts(i,6)
      EFInts(i,1) = ThreeI*XX
      EFInts(i,4) = ThreeI*YY
      EFInts(i,6) = RR
    end do
    nullify(EFInts)
  end if

# ifdef _DEBUGPRINT_

  ! Stored as nZeta,iElem,jElem,iComp

  write(u6,*) ' In EFInt la,lb=',la,lb
  nzab = nZeta*kab
  do iElem=1,nTri_Elem1(la)
    do jElem=1,nTri_Elem1(lb)
      ij = (jElem-1)*nTri_Elem1(la)+iElem
      ip = ip1+nZeta*(ij-1)
      do iComp=1,nComp
        write(Label,'(A,I2,A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',',iComp,') '
        call RecPrt(Label,' ',Array(ip),nZeta,1)
        ip = ip+nzab
      end do
    end do
  end do
# endif

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ip1),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

end subroutine EFInt
