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
! Copyright (C) 1990,1992,1994, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine Tcrtnc_h(Coef1,n1,m1,Coef2,n2,m2,Coef3,n3,m3,Coef4,n4,m4, &
                    ACInt,mabcd,Scrtch,nScr,ACOut, &
                    IndZet,lZeta,IndEta,lEta)
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
!             Modified to back transformation, January '92.            *
!***********************************************************************

use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n1, m1, n2, m2, n3, m3, n4, m4, mabcd, nScr, lZeta, IndZet(lZeta), lEta, IndEta(lEta)
real(kind=wp), intent(in) :: Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4), ACInt(m1*m2*m3*m4,mabcd)
real(kind=wp), intent(out) :: Scrtch(nScr)
! This should be intent(out), but the aliasing/overlap (see above) prevents it
real(kind=wp), intent(_OUT_) :: ACOut(n1*n2*n3*n4,mabcd)
#include "print.fh"
#include "Molcas.fh"
integer(kind=iwp) :: IncVec, ipA2, ipA3, iPrint, iRout, lsize, lZE, nA3, nCache, nVec

iRout = 18
iPrint = nPrint(iRout)
!iPrint = 99

if (iPrint >= 19) call WrCheck('Tcrtnc:P(AB|CD)',ACInt,m1*m2*m3*m4*mabcd)
if (iPrint >= 99) then
  call RecPrt(' In Tcrtnc: P(ab|cd)',' ',ACInt,m1*m2,m3*m4*mabcd)
  call RecPrt(' Coef1',' ',Coef1,n1,m1)
  call RecPrt(' Coef2',' ',Coef2,n2,m2)
  call RecPrt(' Coef3',' ',Coef3,n3,m3)
  call RecPrt(' Coef4',' ',Coef4,n4,m4)
  write(u6,*) n1,n2,n3,n4
end if

! Reduce for contraction matrix
nCache = (3*lCache)/4-n1*m1-n2*m2
lsize = m1*m2+m2*n1
nVec = m3*m4*mabcd
IncVec = min(max(1,nCache/lsize),nVec)
ipA3 = 1
nA3 = nVec*lZeta ! This is the same for the second set!
ipA2 = ipA3+nA3

call TncHlf_h(Coef1,m1,n1,Coef2,m2,n2,lZeta,nVec,IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndZet)

nCache = (3*lCache)/4-n3*m3-n4*m4
lsize = m3*m4+m4*n3
nVec = mabcd*lZeta
IncVec = min(max(1,nCache/lsize),nVec)

lZE = lZeta*lEta
if (mabcd /= 1) then
  call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
  call DGeTMO(ACOut,mabcd,mabcd,lZE,Scrtch,lZE)
  call dcopy_(mabcd*lZE,Scrtch,1,ACOut,1)
else
  call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
end if

if (iPrint >= 59) call RecPrt(' In Tcrtnc: P(ab|cd) ',' ',ACOut,mabcd,lZE)
if (iPrint >= 19) call WrCheck('Tcrtnc:P(ab|cd)',ACOut,lZE*mabcd)

return

end subroutine Tcrtnc_h
