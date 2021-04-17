************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine multdot (wrk,wrksize,
     & nind,mapda,mapia,ssa,mapdb,mapib,ssb,scalar,
     &                     rc)
c
c     this routine do dot product
c     scalar = A(ind)*B(ind)
c
c     nind        - number of indices in both A B (I)
c     mapda        - direct map of A (I)
c     mapia        - inverse map of A (I)
c     ssa        - symmetry state of A (I)
c     mapdb        - direct map of A (I)
c     mapib        - inverse map of A (I)
c     ssb        - symmetry state of B (I)
c     scalar        - final dot product (O)
c     rc        - return (error) code (O)
c
c     N.B. A and B must have same ordering of indexes
c
c     indA     indB     Implemented
c     1        1         Not yet
c     2        2           Yes
c     3        3           Yes
c     4        4           Yes
c     more              No
c
c
#include "wrk.fh"
c
       integer nind,ssa,ssb,rc
       real*8 scalar
c
       integer mapda(0:512,1:6)
       integer mapia(1:8,1:8,1:8)
c
       integer mapdb(0:512,1:6)
       integer mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer symp,symq,symr,iia,iib,possa,possb,nhelp,length
       real*8 scal
c
       rc=0
c
cT    some tests
c
       do 10 nhelp=1,nind
       if (mapda(0,nhelp).ne.mapdb(0,nhelp)) then
c     RC =1 : nonidentical types of indices (NCI/Stup)
       rc=1
       return
       end if
 10     continue
c
       if (mapda(0,5).ne.mapdb(0,5)) then
c     RC =2 : nonidentical number of symmtry blocks in A and B (Stup)
       rc=2
       return
       end if
c
       if (mapda(0,6).ne.mapdb(0,6)) then
c     RC =3 : nonidentical type of A and B (Stup)
       rc=3
       return
       end if
c
       if (ssa.ne.ssb) then
c     RC =4 : nonidentical symmetry state of A and B (Stup)
       rc=4
       return
       end if
c

       if (nind.eq.4) then
c
cI    4 index matrices
c
       scalar=0.0d0
       do 100 iia=1,mapda(0,5)
c
cI.1  def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
       symr=mapda(iia,5)
c     syms is redundant
       possa=mapda(iia,1)
c
cI.2  def parameters of B
       iib=mapib(symp,symq,symr)
       possb=mapdb(iib,1)
c
cI.3  length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),
     &               scal)
       scalar=scalar+scal
       end if
c
 100    continue
c
       else if (nind.eq.3) then
c
cII   3 index matrices
c
       scalar=0.0d0
       do 200 iia=1,mapda(0,5)
c
cII.1 def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
c     symr is redundant
       possa=mapda(iia,1)
c
cII.2 def parameters of B
       iib=mapib(symp,symq,1)
       possb=mapdb(iib,1)
c
cII.3 length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),
     &               scal)
       scalar=scalar+scal
       end if
c
 200    continue
c
       else if (nind.eq.2) then
c
cIII  2 index matrices
c
       scalar=0.0d0
       do 300 iia=1,mapda(0,5)
c
cIII.1def parameters of A
       symp=mapda(iia,3)
c     symq is redundant
       possa=mapda(iia,1)
c
cIII.2def parameters of B
       iib=mapib(symp,1,1)
       possb=mapdb(iib,1)
c
cIII.3length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call mr0u3wt (length,length,length,1,1,wrk(possa),wrk(possb),
     &               scal)
       scalar=scalar+scal
       end if
c
 300    continue
c
       else if (nind.eq.1) then
c
cIV   1 index matrices
c     RC=5 : 1 index matrices (NCI)
       rc=5
       return
c
       else
cV    more than 4 index matrices
c     RC=6 : more than 4 index matrices (NCI)
       rc=6
       return
       end if
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapia)
       end
