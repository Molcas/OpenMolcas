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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      subroutine prettyprint(A,N,M)
* This subroutine pretty prints the NxM matrix A
      implicit none

#include "output.fh"

* Input arguments
      integer N,M
      real*8 A(N,M)

      integer i,j,jStart,jEnd

      do jStart=1,N,5
        jEnd = min(jStart+4, N)
        write(6,'(1x,5i16)')(j,j=jStart,jEnd)
        do i=1,N
          write(6,'(1x,i3,2x,5f16.8)')i,(A(i,j),j=jStart,jEnd)
        end do
        write(6,*)
      end do

      return
      end
