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
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _CHECK_
subroutine Cntrct(First,Coef1,n1,m1,Coef2,n2,m2,Coef3,n3,m3,Coef4,n4,m4,ACInt,mabMin,mabMax,mcdMin,mcdMax,Scrtch,nScrtch,ACOut, &
                  IndZet,lZeta,IndEta,lEta,nComp)
!***********************************************************************
!                                                                      *
! Object: to transform the integrals from primitives to contracted     *
!         basis functions. The subroutine will do both complete and    *
!         incomplete transformations.                                  *
!                                                                      *
! Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!***********************************************************************

use Molcas, only: lCache
use Definitions, only: wp, iwp
#if defined (_DEBUGPRINT_) || defined (_CHECK_)
use Definitions, only: u6
#endif

implicit none
logical(kind=iwp), intent(inout) :: First
integer(kind=iwp), intent(in) :: n1, m1, n2, m2, n3, m3, n4, m4, mabMin, mabMax, mcdMin, mcdMax, nScrtch, lZeta, lEta, nComp, &
                                 IndZet(lZeta), IndEta(lEta)
real(kind=wp), intent(in) :: Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4), &
                             ACInt(n1*n2*n3*n4,nComp*(mabMax-mabMin+1)*(mcdMax-mcdMin+1))
real(kind=wp), intent(inout) :: ACOut(nComp*(mabMax-mabMin+1)*(mcdMax-mcdMin+1),m1*m2*m3*m4)
real(kind=wp), intent(out) :: Scrtch(nScrtch)
integer(kind=iwp) :: IncVec, ipA2, ipA3, lSize, mabcd, ncache_, nVec

mabcd = nComp*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
#ifdef _DEBUGPRINT_
call RecPrt('Cntrct: Coef1',' ',Coef1,n1,m1)
call RecPrt('Cntrct: Coef2',' ',Coef2,n2,m2)
call RecPrt('Cntrct: Coef3',' ',Coef3,n3,m3)
call RecPrt('Cntrct: Coef4',' ',Coef4,n4,m4)
call RecPrt('Cntrct: [a0|c0]',' ',ACInt,lZeta,lEta*mabcd)
write(u6,*) 'IndZet=',IndZet
write(u6,*) 'IndEta=',IndEta
if (.not. First) call RecPrt(' In Cntrct: Partial (a0|c0)',' ',ACOut,mabcd,m1*m2*m3*m4)
#endif
! The idea here is to make the transformation in subblocks
! (size=IncVec) to minimize cache faults. We split the range of
! the compound index such that the contraction coefficients
! (Coef1{n1,m1} & Coef2{n2,m2}), the first quater transformed block
! {n2,IncVec}, for a fixed m1 index, and the half transformed block
! {IncVec,m1,m2} fit into a fixed cache size.

! Reduce for contraction matrices and 3/4th
nCache_ = (3*lCache)/4-n1*m1-n2*m2
! Compute the size of the first quarter and half transformed block.
lsize = n1*n2+n2*m1
! The length of the compound index
nVec = lEta*mabcd
! Compute the size of the increment of which we will run the
! compound index. It is the largest of 1 or the
IncVec = min(max(1,nCache_/lsize),nVec)
! Pointer to the full block of half transformed integrals
! {nVec*m1*m2}
ipA3 = 1
! Pointer to the first quater transformed block {n1*IncVec}
ipA2 = ipA3+nVec*m1*m2
!#define _CHECK_
#ifdef _CHECK_
if (nVec*m1*m2+n2*IncVec > nScrtch) then
  write(u6,*) 'Cntrct: Memory failure 1'
  write(u6,*) 'nVec*m1*m2(A3)=',nVec*m1*m2
  write(u6,*) 'n2*IncVec(A2)=',n2*IncVec
  write(u6,*) 'n2,IncVec=',n2,IncVec
  write(u6,*) 'nVec,lsize=',nVec,lSize
  write(u6,*) 'n1,n2,m1  =',n1,n2,m1
  write(u6,*) 'nScrtch=',nScrtch
  call Abend()
end if
#endif

call CntHlf(Coef1,m1,n1,Coef2,m2,n2,lZeta,nVec,.true.,IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndZet)

#ifdef _DEBUGPRINT_
call RecPrt('Halftransformed',' ',Scrtch(ipA3),nVec,m1*m2)
#endif

nCache_ = (3*lCache)/4-n3*m3-n4*m4
lsize = n3*n4+n4*m3
nVec = mabcd*m1*m2
IncVec = min(max(1,nCache_/lsize),nVec)
#ifdef _CHECK_
if (nVec*m3*m4+n4*IncVec > nScrtch) then
  write(u6,*) 'Cntrct: Memory failure 2'
  call Abend()
end if
#endif

call CntHlf(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,First,IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
First = .false.

#ifdef _DEBUGPRINT_
call RecPrt(' In Cntrct: (a0|c0) ',' ',ACOut,mabcd,m1*m2*m3*m4)
#endif

return

end subroutine Cntrct
