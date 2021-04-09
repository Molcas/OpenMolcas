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
        subroutine transp (AA,BB, dim1,dim2)
c
c this routine do :
c
c AA(a,b) => BB(b,a)
c
        implicit none
        integer dim1,dim2,i,j
        real*8 AA(1:dim1,1:dim2),BB(1:dim2,1:dim1)
c
        do i=1,dim1
        do j=1,dim2
        BB(j,i)=AA(i,j)
        end do
        end do
c
        return
        end
