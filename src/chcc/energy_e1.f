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
        subroutine Energy_E1 (T1n,Fvo,no,nv,E1)
c
c        this routine do:
c        E1 = sum(a,i) T1n(a,i) . Fvo(a,i)
c
c        calculate E1 component of energy
c
        implicit none
        integer no,nv
        real*8 T1n(1)
        real*8 Fvo(1)
        real*8 e1
c
c        help variables
        integer dim1
c
        e1=0.0d0
        dim1=nv*no
        call mr0u3wt (dim1,dim1,dim1,1,1,T1n(1),Fvo(1),e1)
c
        return
        end
