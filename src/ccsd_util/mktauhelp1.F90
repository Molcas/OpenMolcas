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
       subroutine mktauhelp1 (t2,t11,t12,dima,dimb,dimi,dimj,shifta,    &
     &                        shiftb,fact)
!
!     this routine do
!     t2(a,b,i,j) = t2(a,b,i,j) + fact [T11(a,i).T12(b,j)]
!
!     t2        - T2 matrix (I/O)
!     t11     - T1 amplitudes corresponding to spin ia (I)
!     t12     - T1 amplitudes corresponding to spin jb (I)
!     dima    - 1 dimension of T2 (I)
!     dimb    - 2 dimension of T2 (I)
!     dimi    - 3 domension of T2 (I)
!     dimj    - 4 domension of T2 (I)
!     shifta        - number off occ orbitals in spin and symmetry of a (I)
!     shiftb        - number off occ orbitals in spin and symmetry of b (I)
!     fact    - numerical factor (I)
!
!     N.B. symi must be syma and symj must be symb
!
!
       integer dima,dimb,dimi,dimj,shifta,shiftb
       real*8 fact
       real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
       real*8 t11(1:dima,1:dimi)
       real*8 t12(1:dimb,1:dimj)
!
!     help variables
!
       integer i,j,a,b
!
       do 100 j=1,dimj
       do 101 i=1,dimi
       do 102 b=1,dimb
       do 103 a=1,dima
       t2(a,b,i,j)=t2(a,b,i,j)+fact*(t11(a,i)*t12(b,j))
 103    continue
 102    continue
 101    continue
 100    continue
!
       return
! Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(shifta)
         call Unused_integer(shiftb)
       end if
       end
