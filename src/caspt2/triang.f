************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      subroutine triang(nrow,a)

      IMPLICIT NONE

      integer nrow
      real*8 a(nrow**2)

      integer i,j,ij,ji
      real*8 symm

c Convert a square matrix to triangular in-place.

      do 10 i=2,nrow
        do 11 j=1,i-1
          ij=i+(j-1)*nrow
          ji=j+(i-1)*nrow
          symm=0.5D0*(a(ij)+a(ji))
          a(ji)=symm
          a(ij)=symm
  11    continue
  10  continue
      ij=0
      do i=1,nrow
       do j=1,i
        ij=ij+1
        a(ij)=a(j+nrow*(i-1))
       end do
      end do
      return
      end
