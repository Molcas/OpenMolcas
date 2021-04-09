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
       subroutine mv0v1u (length,X,incX,Y,incY,f)
C
C      Y = Y + f . X
C
       implicit none
#include "chcc1.fh"
c
       integer length,incx,incy
       real*8  f
       real*8  X(*)
       real*8  Y(*)
C
C      help variables
C
       integer i
C
       if (mhkey.eq.1) then
c      ESSL
       call daxpy_(length,f,X(1),incx,Y(1),incY)
C
       else
c      Fotran
C
       if (incx*incy.eq.1) then
c
         do 10 i=1, length
         y(i) = y(i)+f*x(i)
 10      continue
       else
         do 20 i=0, length-1
         y(1+i*incy) = y(1+i*incy)+f*x(1+i*incx)
 20      continue
       end if
C
       end if
c
       return
       end
