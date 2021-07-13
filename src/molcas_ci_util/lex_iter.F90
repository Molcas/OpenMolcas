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

subroutine LEX_ITER(N,K,LEX,LAST)

implicit real*8(A-H,O-Z)
dimension LEX(K)
logical LAST

I = K
! Get the first position to be updated
do while ((I > 0) .and. (LEX(I) == N-K+I))
  I = I-1
end do
! If still remaining combinations, update and
! reset all higher positions to lexicographic order
if (I > 0) then
  LEX(I) = LEX(I)+1
  do J=1,K-I
    LEX(I+J) = LEX(I)+J
  end do
! Else, quit finding combinations
else
  LAST = .true.
end if

end subroutine LEX_ITER
