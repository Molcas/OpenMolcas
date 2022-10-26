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

subroutine SQUARE(A,B,ICB,IRB,NROW)

implicit real*8(A-H,O-Z)
dimension A(*), B(*)

! PAM Sep 06: The two special cases here account for
! almost all calls of this code, and written such as
! not to store with non-zero stride--at least some
! rudimentary consideration for cache.
if (ICB == 1) goto 100
if (IRB == 1) goto 200
! General and inefficient code:
IND = 0
do IROW=0,NROW-1
  do ICOL=0,IROW
    IND = IND+1
    B(1+IROW*ICB+ICOL*IRB) = A(IND)
    B(1+ICOL*ICB+IROW*IRB) = A(IND)
  end do
end do
goto 900

100 continue
do IC=0,NROW-1
  do IR=0,IC
    B(1+IR+IC*IRB) = A(1+IR+(IC*(IC+1))/2)
  end do
end do
do IC=0,NROW-2
  do IR=IC+1,NROW-1
    B(1+IR+IC*IRB) = B(1+IC+IR*IRB)
  end do
end do
goto 900

200 continue
do IC=0,NROW-1
  do IR=0,IC
    B(1+IR+IC*ICB) = A(1+IR+(IC*(IC+1))/2)
  end do
end do
do IC=0,NROW-2
  do IR=IC+1,NROW-1
    B(1+IR+IC*ICB) = B(1+IC+IR*ICB)
  end do
end do

900 continue

return

end subroutine SQUARE
