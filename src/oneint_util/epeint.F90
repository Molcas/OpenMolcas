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

subroutine EPEInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, February '91                            *
!***********************************************************************

implicit real*8(A-H,O-Z)
external TNAI, Fake, Cff2D, XRys2D
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
real*8 TC(3), Coori(3,4), Coora(3,4), CoorAC(3,2)
integer iAnga(4), iDCRT(0:7), iStabO(0:7)
logical EQ, NoSpecial
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(2*3,Coora(1,1),1,Coori(1,1),1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1

! Compute FLOP's and size of the work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on.

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CCoor,TC)
  call dcopy_(3,TC,1,CoorAC(1,2),1)
  call dcopy_(3,TC,1,Coori(1,3),1)
  call dcopy_(3,TC,1,Coori(1,4),1)
  call dcopy_(3,TC,1,Coora(1,3),1)
  call dcopy_(3,TC,1,Coora(1,4),1)

  ! Compute primitive integrals before the application of HRR.

  nT = nZeta
  NoSpecial = .true.
  call Rys(iAnga,nt,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coora,CoorAC,mabmin,mabmax,0,0,Array,nArr*nZeta, &
           TNAI,Fake,Cff2D,xRys2D,NoSpecial)

  ! Use the HRR to compute the required primitive integrals.

  call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

  ! Accumulate contributions to the symmetry adapted operator

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipIn),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_integer(nHer)
  call Unused_integer(nOrdOp)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine EPEInt
