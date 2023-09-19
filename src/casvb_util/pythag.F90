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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

real*8 function PYTHAG(A,B)
! FINDS sqrt(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW

real*8 A, B
real*8 P, R, S, T, U

P = max(abs(A),abs(B))
if (P /= 0.0d0) then
  R = (min(abs(A),abs(B))/P)**2
  do
    T = 4.0d0+R
    if (T == 4.0d0) exit
    S = R/T
    U = 1.0d0+2.0d0*S
    P = U*P
    R = (S/U)**2*R
  end do
end if
PYTHAG = P

return

end function PYTHAG
