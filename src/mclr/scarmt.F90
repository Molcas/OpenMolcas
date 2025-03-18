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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine SCARMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,ISCA,SSCA)
! Scatter-add  rows of MATIN to transposed matrix MATUT
!
!  MATUT(J,ISCA(I)) = MATUT(J,ISCA(I)) + SSCA(I)*MATIN(I,J)
!  (if INDEX(I) /= 0)
!
! Jeppe Olsen, Getting LUCIA in shape, Feb1994

implicit real*8(A-H,O-Z)
real*8 MATIN, MATUT
! Input
dimension ISCA(*), SSCA(*), MATIN(NROWIN,NCOLIN)
! Input and Output
dimension MATUT(NCOLUT,NROWUT) ! (MATUT transposed!)

! To get rid of annoying and incorrect compiler warnings
ICINOF = 0

!LBLK = 100
LBLK = 40
NBLK = NCOLIN/LBLK
if (LBLK*NBLK < NCOLIN) NBLK = NBLK+1
do ICINBL=1,NBLK
  if (ICINBL == 1) then
    ICINOF = 1
  else
    ICINOF = ICINOF+LBLK
  end if
  ICINEN = min(ICINOF+LBLK-1,NCOLIN)
  do I=1,NROWIN
    if (ISCA(I) /= 0) then
      S = SSCA(I)
      IROW = ISCA(I)
      do ICOL=ICINOF,ICINEN
        MATUT(ICOL,IROW) = MATUT(ICOL,IROW)+S*MATIN(I,ICOL)
      end do
    end if
  end do
end do

return

end subroutine SCARMT
