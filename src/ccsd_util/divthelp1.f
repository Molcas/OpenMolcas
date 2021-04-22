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
       subroutine divthelp1 (t1,dima,dimi,dp)
c
c     this routine do
c     t1(a,i) = t1(a,i)/(dp(i)-dp(a))
c
c     t1        - T1 matrix (I/O)
c     dima    - v dimension of T1 (I)
c     dimi    - o domension of T1 (I)
c     dp      - diagonal part of Fok (I)
c
c     N.B. Since for T1 i and a are of the same spin, there is no reason
c     to specify spin of dp. It must be automatically the same spin ad i and a.
c
       integer dima,dimi
       real*8 t1(1:dima,1:dimi)
       real*8 dp(*)
c
c     help variables
c
       integer a,i
       real*8 dpi,den
c
       do 100 i=1,dimi
       dpi=dp(i)
       do 101 a=1,dima
c     t1(a,i)=t1(a,i)/(dpi-dp(dimi+a))
c
       den=dpi-dp(dimi+a)
       if (abs(den).lt.1.0d-7) then
       if (abs(t1(a,i)).gt.1.0d-10) then
       t1(a,i)=t1(a,i)/den
       end if
       else
       t1(a,i)=t1(a,i)/den
       end if
c
 101    continue
 100    continue
c
       return
       end
