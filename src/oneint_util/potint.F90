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
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: i, iAnga(4), iDCRT(0:7), iGrid, iPh(3), ipIn, iStabO(0:7), lDCRT, llOper, LmbdT, mabMax, mabMin, nDCRT, &
                     nFLOP, nMem, nOp, nStabO, nT
real(kind=wp) :: Chrg, Coora(3,4), CoorAC(3,2), Coori(3,4), FACT, TC(3)
logical(kind=iwp) :: NoSpecial
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
external :: Fake, TNAI, XCff2D, XRys2D
! Statement function for Cartesian index
integer(kind=iwp) :: nElem, nabSz, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(nOrdOp)

call fzero(rFinal,nZeta*nElem(la)*nElem(lb)*nIC)

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(2*3,Coora,1,Coori,1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1

! Compute FLOP's and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
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

  do i=1,3
    iph(i) = iPhase(i,iDCRT(lDCRT))
  end do
  nOp = NrOpr(iDCRT(lDCRT))

  do iGrid=1,nGrid
    if (iAddPot /= 0) Chrg = ptchrg(iGrid)
    if (Chrg == Zero) cycle

    do i=1,3
      TC(i) = real(iPh(i),kind=wp)*CCoor(i,iGrid)
      CoorAC(i,2) = TC(i)
      Coori(i,3) = TC(i)
      Coori(i,4) = TC(i)
      Coora(i,3) = TC(i)
      Coora(i,4) = TC(i)
    end do

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
