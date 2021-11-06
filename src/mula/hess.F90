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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine ShiftHess(Hess,shift,nDim,nDim2)
!  Purpose:
!    Shifts Hessian to make it positive definite.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

real*8 Hess(nDim,nDim2)
!real*8 U(nDim,nDim2)
!real*8 Hess_lowT(nDim*(nDim+1)/2)
real*8 epsilon
real*8 eigen_min
logical shift
#include "WrkSpc.fh"

call GetMem('U','Allo','Real',ipU,nDim*nDim2)
call GetMem('Hess_low','Allo','Real',ipHess_lowT,nDim*(nDim+1)/2)

! Initialize.
nDim1 = nDim+1
nDimSqr = nDim**2

k = 0
do i=1,nDim
  do j=1,i
    k = k+1
    Work(ipHess_lowT+k-1) = Hess(i,j)
  end do
end do
call dcopy_(nDimSqr,[0.0d0],0,Work(ipU),1)
call dcopy_(nDim,[1.0d0],0,Work(ipU),nDim1)
call Jacob(Work(ipHess_lowT),Work(ipU),nDim,nDim)
call Jacord(Work(ipHess_lowT),Work(ipU),nDim,nDim)
eigen_min = Work(ipHess_lowT)
shift = eigen_min < 0.0d0
if (shift) then
  epsilon = 2*eigen_min
  do i=1,nDim
    Hess(i,i) = Hess(i,i)-epsilon
  end do
end if
call GetMem('U','Free','Real',ipU,nDim*nDim2)
call GetMem('Hess_low','Free','Real',ipHess_lowT,nDim*(nDim+1)/2)

end subroutine ShiftHess
