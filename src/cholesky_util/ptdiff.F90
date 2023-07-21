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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine PtDiff(IDim,Coeff,T,A)
!-----------------------------------------------------------------------
! Function : Calculate derivative
!            j = odd  ; A(I,J) = EXP(-X*Coeff(J+1))
!            j = even ; A(I,J) = -X*Coeff(J-1)*EXP(-X*Coeff(J))
!-----------------------------------------------------------------------

implicit real*8(A-H,O-Z)
real*8 Coeff(40), T(40), A(40,40)

do I=1,IDim
  do J=1,IDim
    X = T(IDim+1-I)
    if (mod(J,2) == 1) then
      A(I,J) = exp(-X*Coeff(J+1))
    else
      A(I,J) = -X*Coeff(J-1)*exp(-X*Coeff(J))
    end if
  end do
end do

return

end subroutine PtDiff
