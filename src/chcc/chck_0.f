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
        subroutine Chck_0 (dim,A)
c
c        check zero
c
        implicit none
        integer dim
        real*8 A(1:dim)
c
c        help var
        integer i,bad
c
        bad=0
        do i=1,dim
          if (abs(A(i)).gt.1.0d-10) then
          bad=bad+1
          end if
        end do
c
        write (6,*) ' Nonzero elements ',bad,dim
c
c
        return
        end
