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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Exp_1(Vector,n1,n2,Array,Fact)
!***********************************************************************
!                                                                      *
! Object: expand an array.                                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 Vector(n1,n2), Array(n1)

do i2=1,n2
  do i1=1,n1
    Vector(i1,i2) = Array(i1)*Fact
  end do
end do

return

end subroutine Exp_1
