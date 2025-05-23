
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Cntrct_mck(First,Coef1,n1,m1,Coef2,n2,m2,Coef3,n3,m3,Coef4,n4,m4,g1In,nGr,Array,nArr,xpre,g1Out,ngr1,nt,IndZet,nZeta, &
                      lZeta,IndEta,nEta,lEta)
!***********************************************************************
!                                                                      *
! Object: to transform the integrals from primitives to contracted     *
!         basis functions. The subroutine will do both complete and    *
!         incomplete transformations.                                  *
!                                                                      *
! Author:      Roland Lindh, Dept. of Theoretical Chemistry, University*
!              of Lund, SWEDEN.                                        *
!                                                                      *
! Modified by: Anders Bernhardsson for direct implementation of the    *
!              calculation of first order derivatives needed for       *
!              response calculation.                                   *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(inout) :: First
integer(kind=iwp), intent(in) :: n1, m1, n2, m2, n3, m3, n4, m4, nGr, nArr, ngr1, nt, nZeta, IndZet(nZeta), lZeta, nEta, &
                                 IndEta(nEta), lEta
real(kind=wp), intent(in) :: Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4), xpre(nt)
real(kind=wp), intent(inout) :: g1In(nT,nGr), Array(nArr), g1Out(nGr1)
#include "Molcas.fh"
integer(kind=iwp) :: iabcdg, IncVec, ip, ipA2, ipA3, lsize, nCache, nVec

!iRout = 18
!iPrint = nPrint(iRout)

!if (iPrint >= 99) call RecPrt(' In Cntrct: ',' ',G1In,nt,nGr)
!if ((iPrint >= 59) .and. (.not. First)) call RecPrt(' In Cntrct: Partial (a0|c0)',' ',G1Out,nGr,m1*m2*m3*m4)

! Cache size is 32 k word (real)

do iabcdg=1,ngr
  G1In(:,iabcdg) = G1In(:,iabcdg)*xpre(:)
end do

! Reduce for contraction matrix
nCache = (3*lCache)/4-n1*m1-n2*m2
lsize = n1*n2+n2*m1
nVec = lEta*nGr
IncVec = min(max(1,nCache/lsize),nVec)
ipA3 = 1
ipA2 = ipA3+nVec*m1*m2
ip = ipA2+n2*IncVec*m1
if (ip > nArr) call Abend()

call CntHlf_mck(Coef1,m1,n1,Coef2,m2,n2,nZeta,lZeta,nVec,.true.,IncVec,G1In,Array(ipA2),Array(ipA3),IndZet)

nCache = (3*lCache)/4-n3*m3-n4*m4
lsize = n3*n4+n4*m3
nVec = nGr*m1*m2
IncVec = min(max(1,nCache/lsize),nVec)
ip = ipA2+n4*IncVec*m3
if (ip > nArr) call Abend()

call CntHlf_mck(Coef3,m3,n3,Coef4,m4,n4,nEta,lEta,nVec,First,IncVec,Array(ipA3),Array(ipA2),G1Out,IndEta)
First = .false.

!if (iPrint >= 59) call RecPrt(' In Cntrct:  ',' ',ACOut,labcdG,m1*m2*m3*m4)

return

end subroutine Cntrct_mck
