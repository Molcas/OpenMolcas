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
       subroutine divthelp2 (t2,dima,dimb,dimi,dimj,dpa,dpb,dpi,dpj,
     &                       shifta,shiftb)
c
c     this routine do
c     t2(a,b,i,j) = t2(a,b,i,j)/(dpi(i)+dpj(j)-dpa(a)-dpb(b))
c
c     t2        - T2 matrix (I/O)
c     dima    - 1 dimension of T2 (I)
c     dimb    - 2 dimension of T2 (I)
c     dimi    - 3 domension of T2 (I)
c     dimj    - 4 domension of T2 (I)
c     dpa     - diagonal part of Fok corresponding to irrep of a (I)
c     dpb     - diagonal part of Fok corresponding to irrep of b (I)
c     dpi     - diagonal part of Fok corresponding to irrep of i (I)
c     dpj     - diagonal part of Fok corresponding to irrep of j (I)
c     shifta        - number off occ orbitals in spin and symmetry of a (I)
c     shiftb        - number off occ orbitals in spin and symmetry of b (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 dpa(*)
       real*8 dpb(*)
       real*8 dpi(*)
       real*8 dpj(*)
c
c     help variables
c
       integer i,j,a,b
       real*8 den,denj,denij,denijb
c
       do 100 j=1,dimj
       denj=dpj(j)
       do 101 i=1,dimi
       denij=denj+dpi(i)
       do 102 b=1,dimb
       denijb=denij-dpb(shiftb+b)
       do 103 a=1,dima
c     t2(a,b,i,j)=t2(a,b,i,j)/(denijb-dpa(shifta+a))
c
       den=denijb-dpa(shifta+a)
       if (abs(den).lt.1.0d-7) then
       if (abs(t2(a,b,i,j)).gt.1.0d-10) then
       t2(a,b,i,j)=t2(a,b,i,j)/den
       end if
       else
       t2(a,b,i,j)=t2(a,b,i,j)/den
       end if
c
 103    continue
 102    continue
 101    continue
 100    continue
c
       return
       end
