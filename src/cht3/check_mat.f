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
        subroutine check_mat(mat,dima,dimb)
c
        implicit none
        integer dima,dimb,i,j
        real*8 mat(dima,dimb)
c
        do i=1,dima
        do j=1,dimb
        if (abs(mat(i,j)).gt.10000) then
          write (6,*) 'i,j,mat(i,j) ',i,j,mat(i,j)
        end if
        end do
        end do
c
        return
        end
