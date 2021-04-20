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
        subroutine GetTauHlp1 (Tau,T1,dima,dimb,adda,addb,no,nv)
c
c        Make Tau for aGrp.ne.bGrp
c
        implicit none
        integer dima,dimb,adda,addb,no,nv
        real*8 Tau(1:dima,1:dimb,1:no,1:no)
        real*8 T1(1:nv,1:no)
c
c       help variables
        integer i,j,a,b
        real*8 c
c
c
        do j=1,no
          do i=1,no
            do b=1,dimb
              c=t1(addb+b,j)
              do a=1,dima
                Tau(a,b,i,j)=Tau(a,b,i,j)+c*t1(adda+a,i)
              end do
            end do
          end do
        end do
c
c
        return
        end
