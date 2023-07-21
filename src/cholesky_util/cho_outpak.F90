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
! Copyright (C) 1973, Nelson H. F. Beebe                               *
!***********************************************************************

subroutine CHO_OUTPAK(AMATRX,NROW,NCTL,LUPRI)
!.......................................................................
!
! OUTPAK PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
!
! TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
!
! ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
!
!        AMATRX(*)...........PACKED MATRIX
!
!        NROW................NUMBER OF ROWS TO BE OUTPUT
!
!        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
!                                                   2 FOR DOUBLE SPACE,
!                                                   3 FOR TRIPLE SPACE.
!
! THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:
!
!        1
!        2    3
!        4    5    6
!        7    8    9   10
!       11   12   13   14   15
!       16   17   18   19   20   21
!       22   23   24   25   26   27   28
!       AND SO ON.
!
! OUTPAK IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
! FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
! COLUMNS.
!
! AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
!..........VERSION = 09/05/73/03
!.......................................................................

implicit real*8(a-h,o-z)
real*8 AMATRX(*)
integer BEGIN
character*1 ASA(3), BLANK, CTL
character PFMT*20, COLUMN*8
parameter(ZERO=0.d00,KCOLP=4,KCOLN=6)
parameter(FFMIN=1.D-3,FFMAX=1.d3)
data COLUMN/'Column  '/,ASA/' ','0','-'/,BLANK/' '/

if (NCTL < 0) then
  KCOL = KCOLN
else
  KCOL = KCOLP
end if
MCTL = abs(NCTL)
if ((MCTL <= 3) .and. (MCTL > 0)) then
  CTL = ASA(MCTL)
else
  CTL = BLANK
end if

J = NROW*(NROW+1)/2
AMAX = ZERO
do I=1,J
  AMAX = max(AMAX,abs(AMATRX(I)))
end do
if (AMAX == ZERO) then
  write(LUPRI,'(/T6,A)') 'Zero matrix.'
  GO TO 200
end if
if ((FFMIN <= AMAX) .and. (AMAX <= FFMAX)) then
  ! use F output format
  PFMT = '(A1,I7,2X,8F15.8)'
else
  ! use 1PD output format
  PFMT = '(A1,I7,2X,1P,8D15.6)'
end if

! LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED

LAST = min(NROW,KCOL)

! BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.

!.....BEGIN NON STANDARD DO LOOP.
BEGIN = 1
1050 NCOL = 1
write(LUPRI,1000) (COLUMN,I,I=BEGIN,LAST)
do K=BEGIN,NROW
  KTOTAL = (K*(K-1))/2+BEGIN-1
  do I=1,NCOL
    if (AMATRX(KTOTAL+I) /= ZERO) GO TO 20
  end do
  GO TO 30
20 write(LUPRI,PFMT) CTL,K,(AMATRX(J+KTOTAL),J=1,NCOL)
30 if (K < (BEGIN+KCOL-1)) NCOL = NCOL+1
end do
LAST = min(LAST+KCOL,NROW)
BEGIN = BEGIN+NCOL
if (BEGIN <= NROW) GO TO 1050
200 continue

return

1000 format(/12X,6(3X,A6,I4,2X),(3X,A6,I4))
!2000 format(A1,'Row',I4,2X,1P,8D15.6)
!2000 format(A1,I7,2X,1P,8D15.6)

end subroutine CHO_OUTPAK
