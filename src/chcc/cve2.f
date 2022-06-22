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
        subroutine CVE2 (V,oe,dima,dimb,adda,addb,no,e2,e2os)
c
c       this routine do:
c       e2  =e2   + sum(a',b',i,j)[(a'i|b'j).{2(a'i|b'j)-(a'j|b'i)}]/Dija'b'
c       e2os=e2os + sum(a',b',i,j)[(a'i|b'j).{ (a'i|b'j)          }]/Dija'b'
c       N.B. Cierny Vypocet E2
c       N.B.2 qvajt odflaknute
c
        implicit none
        integer dima,dimb,adda,addb,no
        real*8 V(1:dima,1:no,1:dimb,1:no)
        real*8 oe(*)
        real*8 e2,e2os
c
c       help variables
        integer i,j,a,b
        real*8 dijab
c
c
        do j=1,no
          do b=1,dimb
            do i=1,no
              do a=1,dima
                dijab=oe(i)+oe(j)-oe(adda+a)-oe(addb+b)
                e2  =e2   +
     &     V(a,i,b,j)*(2.0d0*V(a,i,b,j)-V(a,j,b,i))/dijab
                e2os=e2os +
     &     V(a,i,b,j)*(      V(a,i,b,j)           )/dijab
              end do
            end do
          end do
        end do
c
c
        return
        end
