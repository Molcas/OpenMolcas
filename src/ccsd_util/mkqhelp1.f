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
       subroutine mkqhelp1 (t2,t11,t12,dima,dimb,dimi,dimj,shifta,
     &                      shiftb,fact)
c
c     this routine do
c     t2(a,b,i,j) = fact . t2(a,b,i,j) [T11(a,i).T12(b,j)]
c
c     t2        - T2 matrix (I/O)
c     t11     - T1 amplitudes corresponding to spin ia (I)
c     t12     - T1 amplitudes corresponding to spin jb (I)
c     dima    - 1 dimension of T2 (I)
c     dimb    - 2 dimension of T2 (I)
c     dimi    - 3 domension of T2 (I)
c     dimj    - 4 domension of T2 (I)
c     shifta        - number off occ orbitals in spin and symmetry of a (I)
c     shiftb        - number off occ orbitals in spin and symmetry of b (I)
c     fact    - numerical factor (I)
c
c     N.B. symi must be syma and symj must be symb
c
c
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 fact
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 t11(1:dima,1:dimi)
       real*8 t12(1:dimb,1:dimj)
c
c     help variables
c
       integer i,j,a,b
c
       do 100 j=1,dimj
       do 101 i=1,dimi
       do 102 b=1,dimb
       do 103 a=1,dima
       t2(a,b,i,j)=fact*t2(a,b,i,j)+(t11(a,i)*t12(b,j))
 103    continue
 102    continue
 101    continue
 100    continue
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(shifta)
         call Unused_integer(shiftb)
       end if
       end
