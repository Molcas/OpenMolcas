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
        subroutine T2d_div (T2,OE,dima,dimb,adda,addb,no,nv)
c
c        this routine do:
c        T2(a',b',i,j) = T2(a',b',i,j/(e(i)+e(j)-e(a)-e(b))
c        for aGrp=beGrp, where only a'>=b',i,j are valid,
c        and completed also cases a'<b',i,j
c
c        divison of T2n amplitides by denominator
c
        implicit none
        integer dima,dimb,adda,addb,no,nv
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 OE(1:no+nv)
c
c        help variables
        integer i,j,a,av,b,bv
        real*8 eija
c
        av=no+adda
        bv=no+addb
c
c1        divison by denominators
c
        do j=1,no
        do i=1,no
          do a=1,dima
          eija=OE(i)+OE(j)-OE(av+a)
          do b=1,a
            T2(a,b,i,j)=T2(a,b,i,j)/(eija-OE(bv+b))
          end do
        end do
        end do
        end do
c
c2        completing upper triangle
c
        do j=1,no
        do i=1,no
          do b=2,dima
          do a=1,b-1
             T2(a,b,i,j)=T2(b,a,j,i)
          end do
        end do
        end do
        end do
c
c
        return
        end
