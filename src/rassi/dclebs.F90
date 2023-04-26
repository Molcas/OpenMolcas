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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

function DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,XM3)
! DCLEBS: real Clebsch-Gordan coefficients
! From a modification of Racah''s formula. Coded: Malmqvist 1998
! Note carefully: The input values XJ1..XM3 are REAL, not integers. Half-integer spins are allowed
! Half-integers are assumed exactly represented

use Definitions, only: wp, iwp
use Constants, only: Zero, One

implicit none
real(kind=wp) :: DCLEBS
real(kind=wp), intent(in) :: XJ1,XJ2,XJ3,XM1,XM2,XM3
integer, parameter :: MAXJ=10, MAXF=3*MAXJ+1
integer(kind=iwp), save :: icall = 0
real(kind=wp), save :: DFACT(0:MAXF)
real(kind=wp) :: DF, den, PRE, PRE2, SUMMA, TERM, XJSUM
integer(kind=iwp) :: i, IA1, IA2, IA3, IB1, IB2, IB3, IX, IX1, IX2, IY, IY0, JSUM

if (icall == 0) then
  icall = icall+1
  DF = One
  DFACT(0) = DF
  do i=1,MAXF
    DF = DBLE(i)*DF
    DFACT(i) = DF
  end do
end if

DCLEBS = Zero

XJSUM = XJ1+XJ2+XJ3
JSUM = NINT(XJSUM)
if (XJSUM /= DBLE(JSUM)) return
if (XM1+XM2 /= XM3) return

IA1 = nint(XJ1+XM1)
if (IA1 < 0) return
IB1 = nint(XJ1-XM1)
if (IB1 < 0) return
IA2 = nint(XJ2+XM2)
if (IA2 < 0) return
IB2 = nint(XJ2-XM2)
if (IB2 < 0) return
IA3 = nint(XJ3-XM3)
if (IA3 < 0) return
IB3 = nint(XJ3+XM3)
if (IB3 < 0) return
if (JSUM-IA1-IB1 < 0) return
if (JSUM-IA2-IB2 < 0) return
if (JSUM-IA3-IB3 < 0) return

PRE2 = real(1+IA3+IB3)*DFACT(JSUM-IA1-IB1) &
       *DFACT(JSUM-IA2-IB2)*DFACT(JSUM-IA3-IB3) &
       *DFACT(IA1)*DFACT(IA2)*DFACT(IA3) &
       *DFACT(IB1)*DFACT(IB2)*DFACT(IB3) &
       /DFACT(JSUM+1)
PRE = sqrt(PRE2)

IY0 = (JSUM-IA3-IB3)
IX1 = (IA2+IB1-JSUM)+IB2
IX2 = (IA2+IB1-JSUM)+IA1
IX = MAX(0,IX1,IX2)
IY = MIN(IY0,IB1,IA2)

SUMMA = Zero
do I=IX,IY
  DEN = DFACT(I)*DFACT(I-IX1)*DFACT(I-IX2)*DFACT(IY0-I)*DFACT(IB1-I)*DFACT(IA2-I)
  TERM = One/DEN
  SUMMA = SUMMA+real((-1)**I)*TERM
end do

DCLEBS = PRE*SUMMA

end function DCLEBS
