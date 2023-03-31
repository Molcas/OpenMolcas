!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine unpckhelp5 (a,b,dimp,dimj,dime,jadd,noj,eadd,noe)
!
!     this routine do:
!     b(j,e) = a(pj,qe)-a(qe,pj) for symp=symq
!
       integer dimp,dime,dimj,eadd,noe,jadd,noj
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimj,1:dime)
!
!     help variables
       integer pj,qe,e
!
       do 100 qe=eadd+1,eadd+noe
       e=qe-eadd
       do 101 pj=jadd+1,jadd+noj
       b(pj-jadd,e)=a(pj,qe)-a(qe,pj)
 101    continue
 100    continue
!
       return
       end
