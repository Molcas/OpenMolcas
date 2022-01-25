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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine CZERO2(XX,N,M,NDIM)

implicit real*8(A-H,O-Z)
dimension XX(NDIM,*)
data ZERO/0.0D+00/

do J=1,N
  do I=1,M
    XX(I,J) = ZERO
  end do
end do

return

end subroutine CZERO2
