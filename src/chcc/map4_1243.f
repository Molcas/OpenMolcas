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
        subroutine Map4_1243 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(1243) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1*d2,1:d3,1:d4)
        real*8 B(1:d1*d2,1:d4,1:d3)
c
c       help variables
        integer i12,i3,i4
c
        do i3=1,d3
        do i4=1,d4
        do i12=1,d1*d2
        b(i12,i4,i3)=a(i12,i3,i4)
        end do
        end do
        end do
c
        return
        end
