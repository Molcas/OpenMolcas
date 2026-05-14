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
! Copyright (C) 1990,1992,1994,1996, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _CHECK_
subroutine Tcrtnc(Coef1,n1,m1,Coef2,n2,m2,Coef3,n3,m3,Coef4,n4,m4,ACInt,mabcd,Scrtch,nScr,ACOut,IndZet,lZeta,IndEta,lEta)
!***********************************************************************
!                                                                      *
! Object: to transform the integrals from primitives to contracted     *
!         basis functions. The subroutine will do both complete and    *
!         incomplete transformations.                                  *
!                                                                      *
!         Observe that ACInt and ACOut may overlap!!!!                 *
!         (which is against Fortran standard, so this should be fixed) *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to back transformation, January '92.            *
!***********************************************************************

use Molcas, only: lCache
use Definitions, only: wp, iwp
#if defined (_DEBUGPRINT_) || defined (_CHECK_)
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n1, m1, n2, m2, n3, m3, n4, m4, mabcd, nScr, lZeta, IndZet(lZeta), lEta, IndEta(lEta)
real(kind=wp), intent(in) :: Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4), ACInt(m1*m2*m3*m4,mabcd)
real(kind=wp), intent(out) :: Scrtch(nScr)
! FIXME: This should be intent(out), but the aliasing/overlap (see above) prevents it
real(kind=wp), intent(_OUT_) :: ACOut(lZeta*lEta,mabcd)
integer(kind=iwp) :: nCache, lsize, nVec, IncVec, ipA2, ipA3

#ifdef _DEBUGPRINT_
call WrCheck('Tcrtnc:P(AB|CD)',ACInt,m1*m2*m3*m4*mabcd)
call RecPrt(' In Tcrtnc: P(ab|cd)',' ',ACInt,mabcd,m1*m2*m3*m4)
call RecPrt(' Coef1',' ',Coef1,n1,m1)
call RecPrt(' Coef2',' ',Coef2,n2,m2)
call RecPrt(' Coef3',' ',Coef3,n3,m3)
call RecPrt(' Coef4',' ',Coef4,n4,m4)
write(u6,*) n1,n2,n3,n4
#endif

! Reduce for contraction matrix
nCache = (3*lCache)/4-n3*m3-n4*m4
lsize = m3*m4+m3*n4
nVec = m1*m2*mabcd
IncVec = min(max(1,nCache/lsize),nVec)
ipA3 = 1
ipA2 = ipA3+nVec*lEta

#ifdef _CHECK_
if (nVec*lEta+n4*m3*IncVec > nScr) then
  write(u6,*) 'Tcrtnc: Memory failure 1'
  write(u6,*) 'n4*IncVec*m3(A2)=',n4*IncVec*m3
  write(u6,*) 'nVec*lEta(A3)=',nVec*lEta
  write(u6,*) 'n4,IndVec,m3=',n4,IndVec,m3
  write(u6,*) 'nVec,lEta=',nVec,lEta
  write(u6,*) 'nScr=',nScr
  call Abend()
end if
#endif

call TncHlf(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndEta)

nCache = (3*lCache)/4-n1*m1-n2*m2
lsize = m1*m2+m1*n2
nVec = mabcd*lEta
IncVec = min(max(1,nCache/lsize),nVec)

#ifdef _CHECK_
if (nVec*m1*m2+n2*IncVec*m1 > nScr) then
  write(u6,*) 'Tcrtnc: Memory failure 2'
  write(u6,*) 'nVec*m1*m2(A1)=',nVec*m1*m2
  write(u6,*) 'n2*IncVec*m1(A2)=',n2*IncVec*m1
  write(u6,*) 'nVec,m1,m2=',nVec,m1,m2
  write(u6,*) 'n2,IncVec,m1=',n2,IncVec,m1
  write(u6,*) 'nScr=',nScr
  call Abend()
end if
#endif

call TncHlf(Coef1,m1,n1,Coef2,m2,n2,lZeta,nVec,IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndZet)

#ifdef _DEBUGPRINT_
call RecPrt(' In Tcrtnc: P(ab|cd) ',' ',ACOut,mabcd,lZeta*lEta)
call WrCheck('Tcrtnc:P(ab|cd)',ACOut,lZeta*lEta*mabcd)
#endif

end subroutine Tcrtnc
