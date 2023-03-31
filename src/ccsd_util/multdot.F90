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
       subroutine multdot (wrk,wrksize,                                 &
     & nind,mapda,mapia,ssa,mapdb,mapib,ssb,scalar,                     &
     &                     rc)
!
!     this routine do dot product
!     scalar = A(ind)*B(ind)
!
!     nind        - number of indices in both A B (I)
!     mapda        - direct map of A (I)
!     mapia        - inverse map of A (I)
!     ssa        - symmetry state of A (I)
!     mapdb        - direct map of A (I)
!     mapib        - inverse map of A (I)
!     ssb        - symmetry state of B (I)
!     scalar        - final dot product (O)
!     rc        - return (error) code (O)
!
!     N.B. A and B must have same ordering of indexes
!
!     indA     indB     Implemented
!     1        1         Not yet
!     2        2           Yes
!     3        3           Yes
!     4        4           Yes
!     more              No
!
!
#include "wrk.fh"
!
       integer nind,ssa,ssb,rc
       real*8 scalar
!
       integer mapda(0:512,1:6)
       integer mapia(1:8,1:8,1:8)
!
       integer mapdb(0:512,1:6)
       integer mapib(1:8,1:8,1:8)
!
!     help variables
!
       integer symp,symq,symr,iia,iib,possa,possb,nhelp,length
       real*8 scal
!
       rc=0
!
!T    some tests
!
       do 10 nhelp=1,nind
       if (mapda(0,nhelp).ne.mapdb(0,nhelp)) then
!     RC =1 : nonidentical types of indices (NCI/Stup)
       rc=1
       return
       end if
 10     continue
!
       if (mapda(0,5).ne.mapdb(0,5)) then
!     RC =2 : nonidentical number of symmtry blocks in A and B (Stup)
       rc=2
       return
       end if
!
       if (mapda(0,6).ne.mapdb(0,6)) then
!     RC =3 : nonidentical type of A and B (Stup)
       rc=3
       return
       end if
!
       if (ssa.ne.ssb) then
!     RC =4 : nonidentical symmetry state of A and B (Stup)
       rc=4
       return
       end if
!

       if (nind.eq.4) then
!
!I    4 index matrices
!
       scalar=0.0d0
       do 100 iia=1,mapda(0,5)
!
!I.1  def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
       symr=mapda(iia,5)
!     syms is redundant
       possa=mapda(iia,1)
!
!I.2  def parameters of B
       iib=mapib(symp,symq,symr)
       possb=mapdb(iib,1)
!
!I.3  length must be common for both A and B
       length=mapda(iia,2)
!
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),    &
     &               scal)
       scalar=scalar+scal
       end if
!
 100    continue
!
       else if (nind.eq.3) then
!
!II   3 index matrices
!
       scalar=0.0d0
       do 200 iia=1,mapda(0,5)
!
!II.1 def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
!     symr is redundant
       possa=mapda(iia,1)
!
!II.2 def parameters of B
       iib=mapib(symp,symq,1)
       possb=mapdb(iib,1)
!
!II.3 length must be common for both A and B
       length=mapda(iia,2)
!
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),    &
     &               scal)
       scalar=scalar+scal
       end if
!
 200    continue
!
       else if (nind.eq.2) then
!
!III  2 index matrices
!
       scalar=0.0d0
       do 300 iia=1,mapda(0,5)
!
!III.1def parameters of A
       symp=mapda(iia,3)
!     symq is redundant
       possa=mapda(iia,1)
!
!III.2def parameters of B
       iib=mapib(symp,1,1)
       possb=mapdb(iib,1)
!
!III.3length must be common for both A and B
       length=mapda(iia,2)
!
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),    &
     &               scal)
       scalar=scalar+scal
       end if
!
 300    continue
!
       else if (nind.eq.1) then
!
!IV   1 index matrices
!     RC=5 : 1 index matrices (NCI)
       rc=5
       return
!
       else
!V    more than 4 index matrices
!     RC=6 : more than 4 index matrices (NCI)
       rc=6
       return
       end if
!
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapia)
       end
