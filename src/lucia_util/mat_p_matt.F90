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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************

subroutine MAT_P_MATT(A,B,NR,NC,COEF)
! A(I,J) = A(I,J) + Coef*B(J,I)

implicit real*8(A-H,O-Z)
! Input and output
dimension A(NR,NC)
! Input
dimension B(NC,NR)

do J=1,NC
  do I=1,NR
    A(I,J) = A(I,J)+COEF*B(J,I)
  end do
end do

end subroutine MAT_P_MATT
