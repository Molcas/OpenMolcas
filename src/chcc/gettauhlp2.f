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
         subroutine GetTauHlp2 (Tau,T1,dima,adda,no,nv)
c
c        Make Tau for aGrp.eq.bGrp
c
        implicit none
         integer dima,adda,no,nv
        real*8 Tau(1:dima*(dima+1)/2,1:no,1:no)
        real*8 T1(1:nv,1:no)
c
c       help variables
        integer i,j,a,b,ab
        real*8 c
c
c
        do j=1,no
          do i=1,no
            ab=0
            do a=1,dima
              c=t1(adda+a,i)
              do b=1,a
                ab=ab+1
                Tau(ab,i,j)=Tau(ab,i,j)+c*t1(adda+b,j)
              end do
            end do
          end do
        end do
c
c
        return
        end
