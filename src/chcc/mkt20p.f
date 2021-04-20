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
        subroutine MkT20p (T2,V,oe,dima,adda,no)
c
c       this routine do:
c       T2(a'b',i,j)=(a'i|b'j)/Dija'b'
c        for aGrp.eq.bGrp
c       N.B.2 qvajt odflaknute, neurobene
c
        implicit none
        integer dima,adda,no
        real*8 V(1:dima,1:no,1:dima,1:no)
        real*8 T2(1:dima*(dima+1)/2,1:no,1:no)
        real*8 oe(*)
c
c       help variables
        integer i,j,a,b,ab
        real*8 dijab
c
c
        do j=1,no
          do i=1,no
            ab=0
            do a=1,dima
            do b=1,a
              ab=ab+1
              dijab=oe(i)+oe(j)-oe(adda+a)-oe(adda+b)
              T2(ab,i,j)=V(a,i,b,j)/dijab
            end do
            end do
          end do
        end do
c
c
        return
        end
