!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine mv0u (length,X,incX,Y,incY)
!
!      Y = X
!
       implicit none
#include "chcc1.fh"
!
       integer length,incx,incy
       real*8  X(*)
       real*8  Y(*)
!
!      help variables
!
       integer i
!
        if (mhkey.eq.1) then
!      ESSL
       call dcopy_(length,X(1),incx,Y(1),incY)
!
       else
!      Fotran
!
       if (incx*incy.eq.1) then
!
         do 10 i=1, length
         y(i) = x(i)
 10      continue
       else
         do 20 i=0, length-1
         y(1+i*incy) = x(1+i*incx)
 20      continue
       end if
!
       end if
!
       return
       end
