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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************

function EQ(A,B)
!***********************************************************************
!                                                                      *
! Object: to return the value .true. if A and B are the same centers.  *
!                                                                      *
!***********************************************************************

implicit none
logical EQ
real*8, intent(In) :: A(3), B(3)

EQ = (A(1) == B(1)) .and. (A(2) == B(2)) .and. (A(3) == B(3))

return

end function EQ
