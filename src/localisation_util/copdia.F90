!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************
! Version of Oct 21

subroutine COPDIA(A,VEC,NDIM,IPACK)
! Copy diagonal of matrix A into vector VEC
!
!   IPACK = 0 : Full matrix
!   IPACK = 1 : Lower triangular matrix

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: A(*), VEC(*)
integer(kind=iwp) :: NDIM, IPACK
integer(kind=iwp) :: I, ip_CPDIA

#include "WrkSpc.fh"
!VVP
!Workaround for dummy aliasing
call GETMEM('CPDIA','ALLO','REAL',ip_CPDIA,NDIM)
if (IPACK == 0) then
  do I=1,NDIM
    Work(ip_CPDIA+I-1) = A((I-1)*NDIM+I)
  end do
else
  do I=1,NDIM
    Work(ip_CPDIA+I-1) = A(I*(I+1)/2)
  end do
end if

call DCOPY_(NDIM,Work(ip_CPDIA),1,VEC,1)
call GETMEM('CPDIA','FREE','REAL',ip_CPDIA,NDIM)

return

end subroutine COPDIA
