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
! Copyright (C) 1999, Markus P. Fuelscher                              *
!               2012, Giovanni Li Manni                                *
!***********************************************************************

subroutine Sort_Cdet(N,Index,X)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Sort CI-coefficients in determinant basis and change phase       *
!                                                                      *
!     calling arguments:                                               *
!     N       : integer                                                *
!               dimension of the CI-vector                             *
!     Index   : integer                                                *
!               reordering indices                                     *
!     X       : real*8                                                 *
!               CI-vector                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     G. Li Manni on 03 Feb 2012 on a desperate situation for saving   *
!     time!                                                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1999                                 *
!                                                                      *
! It was a great piece of code for those times when memory was a very  *
! big issue! Nowadays only time is important!                          *
!***********************************************************************
!do i=1,N                                                              *
!  i_old = i                                                           *
!  i_new = abs(Index(i_old))                                           *
!  do while (i_new > i)                                                *
!    i_old = i_new                                                     *
!    i_new = abs(Index(i_old))                                         *
!  end do                                                              *
!  if (i_new == i) then                                                *
!    i_old = i                                                         *
!    X_old = X(i_old)                                                  *
!    i_new = abs(Index(i_old))                                         *
!    X_new = X(i_new)                                                  *
!    alpha = dble(sign(1,Index(i_old)))                                *
!    do while (i_new > i)                                              *
!      X(i_new) = alpha*X_old                                          *
!      i_old = i_new                                                   *
!      X_old = X_new                                                   *
!      i_new = abs(Index(i_old))                                       *
!      X_new = X(i_new)                                                *
!      alpha = dble(sign(1,Index(i_old)))                              *
!    end do                                                            *
!    X(i) = alpha*X_old                                                *
!  end if                                                              *
!end do                                                                *
!***********************************************************************

implicit real*8(A-H,O-Z)
dimension index(N), X(N)
intrinsic sign
#include "WrkSpc.fh"

call getmem('ReordSDs','Allo','Real',iReoSDs,N)
do i=1,N
  i_new = abs(index(i))
  alpha = dble(sign(1,index(i)))
  Work(iReoSDs+i_new-1) = alpha*X(i)
end do
call dcopy_(N,Work(iReoSDs),1,X,1)
call getmem('ReordSDs','Free','Real',iReoSDs,N)

return

end subroutine Sort_Cdet
