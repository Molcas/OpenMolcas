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

subroutine MINV(ARRAY,ARRINV,DET,NDIM)

implicit real*8(a-h,o-z)
real*8 ARRAY(NDIM,NDIM), ARRINV(NDIM,NDIM)
#include "WrkSpc.fh"

call Allocate_Work(ipA,NDIM**2)
call Allocate_Work(ipB,NDIM**2)
call Allocate_Work(ipBUF,NDIM)
call Allocate_iWork(IPIV,NDIM)
call Allocate_iWork(JPIV,NDIM)

call MINV_INNER(ARRAY,ARRINV,DET,NDIM,Work(ipA),Work(ipBUF),Work(ipB),iWork(IPIV),iWork(JPIV))

call Free_iWork(JPIV)
call Free_iWork(IPIV)
call Free_Work(ipBUF)
call Free_Work(ipB)
call Free_Work(ipA)

return

end subroutine MINV
