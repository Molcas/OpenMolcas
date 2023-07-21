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
! Copyright (C) 1971, Nelson H. F. Beebe                               *
!***********************************************************************

subroutine CHO_OUTPUT(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI)
!.......................................................................
!
! OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;
!
!        AMATRX(',').........MATRIX TO BE OUTPUT
!
!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
!
!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
!
!        ROWDIM..............ROW DIMENSION OF AMATRX(',')
!
!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
!
!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!
! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.
!
! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971
!
!.......................................................................

implicit real*8(a-h,o-z)
integer ROWLOW, ROWHI, COLLOW, COLHI, ROWDIM, COLDIM, BEGIN, KCOL
real*8 AMATRX(ROWDIM,COLDIM)
character*1 ASA(3), BLANK, CTL
character PFMT*20, COLUMN*8
parameter(ZERO=0.d00,KCOLP=4,KCOLN=6)
parameter(FFMIN=1.D-3,FFMAX=1.d3)
data COLUMN/'Column  '/,BLANK/' '/,ASA/' ','0','-'/

if (ROWHI < ROWLOW) GO TO 3
if (COLHI < COLLOW) GO TO 3

AMAX = ZERO
do J=COLLOW,COLHI
  do I=ROWLOW,ROWHI
    AMAX = max(AMAX,abs(AMATRX(I,J)))
  end do
end do
if (AMAX == ZERO) then
  write(LUPRI,'(/T6,A)') 'Zero matrix.'
  GO TO 3
end if
if ((FFMIN <= AMAX) .and. (AMAX <= FFMAX)) then
  ! use F output format
  PFMT = '(A1,I7,2X,8F15.8)'
else
  ! use 1PD output format
  PFMT = '(A1,I7,2X,1P,8D15.6)'
end if

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

LAST = min(COLHI,COLLOW+KCOL-1)
do BEGIN=COLLOW,COLHI,KCOL
  write(LUPRI,1000) (COLUMN,I,I=BEGIN,LAST)
  do K=ROWLOW,ROWHI
    do I=BEGIN,LAST
      if (AMATRX(K,I) /= ZERO) GO TO 5
    end do
    GO TO 1
5   write(LUPRI,PFMT) CTL,K,(AMATRX(K,I),I=BEGIN,LAST)
1   continue
  end do
  LAST = min(LAST+KCOL,COLHI)
end do

3 return

1000 format(/12X,6(3X,A6,I4,2X),(3X,A6,I4))
!2000 format(A1,'Row',I4,2X,1P,8D15.6)
!2000 format(A1,I7,2X,1P,8D15.6)

end subroutine CHO_OUTPUT
