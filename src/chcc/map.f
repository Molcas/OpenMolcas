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
c
c       this file contains:
c
c       Map4_1342
c       Map4_1423
c       Map4_1243
c       Map4_1432
c       Map4_1324
c
c       Map4_2143
c       Map4_2314
c
c       Map4_3124
c       Map4_3142
c       Map4_3214
c       Map4_3412
c       Map4_3421
c
c       Map4_4123
c       Map4_4312
c
c       Map3_132
c       Map3_213
c       Map3_321
c
c       Map2_21
c
c
c       --------------------------------
c
        subroutine Map4_1342 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(1423) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d1,1:d4,1:d2,1:d3)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i3=1,d3
        do i2=1,d2
        do i4=1,d4
        do i1=1,d1
        b(i1,i4,i2,i3)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_1423 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(1342) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d1,1:d3,1:d4,1:d2)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i2=1,d2
        do i4=1,d4
        do i3=1,d3
        do i1=1,d1
        b(i1,i3,i4,i2)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
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
c
c       --------------------------------
c
        subroutine Map4_1432 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(1432) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d1,1:d4,1:d3,1:d2)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i2=1,d2
        do i3=1,d3
        do i4=1,d4
        do i1=1,d1
        b(i1,i4,i3,i2)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_1324 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(1324) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d1,1:d3,1:d2,1:d4)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i4=1,d4
        do i2=1,d2
        do i3=1,d3
        do i1=1,d1
        b(i1,i3,i2,i4)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_2143 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(2143) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d2,1:d1,1:d4,1:d3)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i3=1,d3
        do i4=1,d4
        do i1=1,d1
        do i2=1,d2
        b(i2,i1,i4,i3)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_2314 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(3124) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d3,1:d1,1:d2,1:d4)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i4=1,d4
        do i2=1,d2
        do i1=1,d1
        do i3=1,d3
        b(i3,i1,i2,i4)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_3124 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(2314) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d2,1:d3,1:d1,1:d4)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i4=1,d4
        do i1=1,d1
        do i3=1,d3
        do i2=1,d2
        b(i2,i3,i1,i4)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_3142 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(2413) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d2,1:d4,1:d1,1:d3)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i3=1,d3
        do i1=1,d1
        do i4=1,d4
        do i2=1,d2
        b(i2,i4,i1,i3)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_3214 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(3214) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d3,1:d2,1:d1,1:d4)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i4=1,d4
        do i1=1,d1
        do i2=1,d2
        do i3=1,d3
        b(i3,i2,i1,i4)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_3412 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(3412) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d3,1:d4,1:d1,1:d2)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i2=1,d2
        do i1=1,d1
        do i4=1,d4
        do i3=1,d3
        b(i3,i4,i1,i2)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_3421 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(4312) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d4,1:d3,1:d1,1:d2)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i2=1,d2
        do i1=1,d1
        do i3=1,d3
        do i4=1,d4
        b(i4,i3,i1,i2)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_4123 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(2341) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d2,1:d3,1:d4,1:d1)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i1=1,d1
        do i4=1,d4
        do i3=1,d3
        do i2=1,d2
        b(i2,i3,i4,i1)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map4_4312 (A,B,d1,d2,d3,d4)
c
c       this routine do:
c       map B(3421) <- A(1234)
c
        implicit none
        integer d1,d2,d3,d4
        real*8 A(1:d1,1:d2,1:d3,1:d4)
        real*8 B(1:d3,1:d4,1:d2,1:d1)
c
c       help variables
        integer i1,i2,i3,i4
c
        do i1=1,d1
        do i2=1,d2
        do i4=1,d4
        do i3=1,d3
        b(i3,i4,i2,i1)=a(i1,i2,i3,i4)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map3_132 (A,B,d1,d2,d3)
c
c       this routine do:
c       map B(132) <- A(123)
c
        implicit none
        integer d1,d2,d3
        real*8 A(1:d1,1:d2,1:d3)
        real*8 B(1:d1,1:d3,1:d2)
c
c       help variables
        integer i1,i2,i3
c
        do i2=1,d2
        do i3=1,d3
        do i1=1,d1
        b(i1,i3,i2)=a(i1,i2,i3)
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map3_213 (A,B,d1,d2,d3)
c
c       this routine do:
c       map B(213) <- A(123)
c
        implicit none
        integer d1,d2,d3
        real*8 A(1:d1,1:d2,1:d3)
        real*8 B(1:d2,1:d1,1:d3)
c
c       help variables
        integer i1,i2,i3
c
        do i3=1,d3
        do i2=1,d2
        do i1=1,d1
        b(i2,i1,i3)=a(i1,i2,i3)
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map3_321 (A,B,d1,d2,d3)
c
c       this routine do:
c       map B(321) <- A(123)
c
        implicit none
        integer d1,d2,d3
        real*8 A(1:d1,1:d2,1:d3)
        real*8 B(1:d3,1:d2,1:d1)
c
c       help variables
        integer i1,i2,i3
c
        do i1=1,d1
        do i2=1,d2
        do i3=1,d3
        b(i3,i2,i1)=a(i1,i2,i3)
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Map2_21 (A,B,d1,d2)
c
c       this routine do:
c       map B(21) <- A(12)
c
        implicit none
        integer d1,d2
        real*8 A(1:d1,1:d2)
        real*8 B(1:d2,1:d1)
c
c       help variables
        integer i1,i2
c
        do i1=1,d1
        do i2=1,d2
        b(i2,i1)=a(i1,i2)
        end do
        end do
c
        return
        end
