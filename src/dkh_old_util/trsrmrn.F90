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
! Copyright (C) 1995, Bernd Artur Hess                                 *
!***********************************************************************

subroutine TrSmrN(A,BA,BB,C,NA,NB,H,W)
! TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION IN B.
! RESULT IS IN C

implicit real*8(A-H,O-Z)
dimension A(NA,NB), BA(NA,NA), BB(NB,NB), C(NA,NB), H(NA,NB), W(NA,NB)

call DZero(C,NA*NB)
call DZero(H,NA*NB)
call dcopy_(NA*NB,A,1,W,1)

do I=1,NA
  do L=1,NB
    do K=1,NA
      H(I,L) = BA(K,I)*W(K,L)+H(I,L)
    end do
  end do
end do

do I=1,NA
  do J=1,NB
    do L=1,NB
      C(I,J) = H(I,L)*BB(L,J)+C(I,J)
    end do
  end do
end do

return

end subroutine TrSmrN
