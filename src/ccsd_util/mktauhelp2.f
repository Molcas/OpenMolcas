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
       subroutine mktauhelp2 (t2,t1,dimab,dimij,dima,dimi,shift,fact)
c
c     this routine do
c     t2(ab,ij) = t2(ab,ij) + fact* (T1(a,i).t1(b,j)-t1(b,i).t1(a,j))
c     for spin and symmetry of all indices equal
c
c     t2        - T2 matrix (I/O)
c     t1      - t1 matrix (I)
c     dimab   - 1 dimension of T2 (I)
c     dimij   - 3 domension of T2 (I)
c     dima    - number of a in this spin and symm of a (I)
c     dimi    - number of a in this spin and symm of i (I)
c     shift         - number off occ orbitals in spin and symmetry of a (I)
c     fact    - numerical factor (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dimab,dimij,dima,dimi,shift
       real*8 fact
       real*8 t2(1:dimab,1:dimij)
       real*8 t1(1:dima,1:dimi)
c
c     help variables
c
       integer i,j,a,b,ij,ab
c
       ij=0
       do 100 i=2,dimi
       do 101 j=1,i-1
       ij=ij+1
c
       ab=0
       do 102 a=2,dima
       do 103 b=1,a-1
       ab=ab+1
       t2(ab,ij)=t2(ab,ij)+fact*(t1(a,i)*t1(b,j)-t1(b,i)*t1(a,j))
c
 103    continue
 102    continue
 101    continue
 100    continue
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(shift)
       end
