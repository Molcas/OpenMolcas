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
        subroutine RNFill (length,A,c)
c
c       fill an array with random numbers in interval (-c,c)
c
        implicit none
        integer length
        real*8 c
        real*8 A(1:length)
c
c       help variables
        integer i
c
c
        do i=1,length
c          A(i)=c*(srand()-0.5d0)
           A(i)=(1.0d-7)*i
        end do
c
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_real(c)
        end
