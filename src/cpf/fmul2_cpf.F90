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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine FMUL2_CPF(A,B,C,NROW,NCOL,N)

implicit real*8(A-H,O-Z)
dimension A(NROW,N), B(NCOL,N), CJ(200)
dimension C(NROW,NCOL)

if (nRow > 200) then
  write(6,*)
  call XFLUSH(6)
  write(6,*) ' *** Error in Subroutine FMUL2_CPF ***'
  call XFLUSH(6)
  write(6,*) ' row dimension exceeds local buffer size'
  call XFLUSH(6)
  write(6,*)
  call XFLUSH(6)
  call Abend()
end if

do J=1,NCOL
  do I=1,NROW
    CJ(I) = 0.0d0
  end do
  if (J == NCOL) GO TO 16
  J1 = J+1
  do K=1,N
    FAC = B(J,K)
    if (FAC == 0.0) GO TO 20
    do I=J1,NROW
      CJ(I) = CJ(I)+FAC*A(I,K)
    end do
20  continue
  end do
16 do I=1,NROW
    C(I,J) = CJ(I)
  end do
end do

return

end subroutine FMUL2_CPF
