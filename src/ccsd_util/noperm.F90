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
       subroutine noperm (wrk,wrksize,                                  &
     & mapda,mapia,mapdb,mapib,poss0,posst)
!
!     realize mapping without permutation
!     define mapd,mapi
!
#include "ccsd1.fh"
#include "wrk.fh"
!
       integer poss0,posst
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
!
!     help variables
!
       integer ib,nhelp,i,j,k
!
!     def mapib
!
       do 10 k=1,nsym
       do 11 j=1,nsym
       do 12 i=1,nsym
       mapib(i,j,k)=mapia(i,j,k)
 12     continue
 11     continue
 10     continue
!
!     def initial values
!
       do nhelp=1,6
       mapdb(0,nhelp)=mapda(0,nhelp)
       end do
!
       posst=poss0
       do 100 ib=1,mapda(0,5)
       do nhelp=2,6
       mapdb(ib,nhelp)=mapda(ib,nhelp)
       end do
       mapdb(ib,1)=posst
       posst=posst+mapdb(ib,2)
!
       call map11 (wrk(mapda(ib,1)),wrk(mapdb(ib,1)),mapda(ib,2),1)
!
 100    continue
!
       return
       end
