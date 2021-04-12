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
        subroutine Energy_E2od (V,Tau,e,eos,dima,dimb,no)
c
c        this routine calc:
c        E2 = sum(a,b,i,j) (2V(ai|bj)-V(aj|bi)) . Tau(a,b,i,j)
c        E2os=sum(a,b,i,j)  V(ai|bj)            . Tau(a,b,i,j)
c        where
c           E2    - complete E2 component of energy
c           E2os  - other spin E2 component of energy
c
        implicit none
        integer dima,dimb,no
        real*8 V(1:dima,1:no,1:dimb,1:no)
        real*8 Tau(1:dima,1:dimb,1:no,1:no)
        real*8 e,eos
c
c        help variables
        integer a,b,i,j
c
        e=0.0d0
        eos=0.0d0
c
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
           e=e+(2.0d0*V(a,i,b,j)-V(a,j,b,i))*Tau(a,b,i,j)
           eos=eos+V(a,i,b,j)*Tau(a,b,i,j)
        end do
        end do
        end do
        end do
c
        return
        end
