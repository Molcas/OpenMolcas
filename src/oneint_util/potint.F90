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

subroutine PotInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of potential integrals    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!***********************************************************************

use Phase_Info, only: iPhase
use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, rys2d_kernel, tval_kernel
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: iAnga(4), iDCRT(0:7), iGrid, iPh(3), ipIn, iStabO(0:7), lDCRT, llOper, LmbdT, mabMax, mabMin, nDCRT, nFLOP, &
                     nMem, nOp, nStabO, nT
real(kind=wp) :: Chrg, Coora(3,4), CoorAC(3,2), Coori(3,4), FACT, TC(3)
logical(kind=iwp) :: NoSpecial
procedure(cff2d_kernel) :: XCff2D
procedure(modu2_kernel) :: Fake
procedure(rys2d_kernel) :: XRys2D
procedure(tval_kernel) :: TNAI
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(nOrdOp)

rFinal(:,:,:,:) = Zero

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
Coora(:,1) = A
Coora(:,2) = RB
Coori(:,1:2) = Coora(:,1:2)
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)

! Compute FLOP's and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if

llOper = lOper(1)

! Loop over grid

! Find the DCR for M and S

call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
!Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)
FACT = One

nT = nZeta
Chrg = -One
NoSpecial = .true.

do lDCRT=0,nDCRT-1

  iPh(:) = iPhase(:,iDCRT(lDCRT))
  nOp = NrOpr(iDCRT(lDCRT))

  do iGrid=1,nGrid
    if (iAddPot /= 0) Chrg = ptchrg(iGrid)
    if (Chrg == Zero) cycle

    TC(:) = real(iPh(:),kind=wp)*CoorO(:,iGrid)
    CoorAC(:,2) = TC(:)
    Coori(:,3) = TC(:)
    Coori(:,4) = TC(:)
    Coora(:,3) = TC(:)
    Coora(:,4) = TC(:)

    ! Compute integrals with the Rys quadrature.

    call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coora,CoorAC,mabmin,mabmax,0,0,Array, &
             nArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

    ! Use the HRR to compute the required primitive integrals.

    call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

    ! Accumulate contributions to the symmetry adapted operator

    call SymAdO(Array(ipIn),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,-Fact*Chrg)
  end do
end do

return

end subroutine PotInt
