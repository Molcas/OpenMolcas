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

subroutine TPQSET(ICASE,TPQ,IP)

use cpf_global, only: INCPF, IPRINT, IRC, IREF0, ISDCI, LN, LWSP, N
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), IP
real(kind=wp), intent(_OUT_) :: TPQ(*)
integer(kind=iwp) :: I, II, II1, IINT, IIOR, IJ, IK, IL, IOCR(100), IQ, JJ, JOJ, NI, NJ
real(kind=wp) :: DIK, DIL, DJK, DJL
integer(kind=iwp), external :: ICUNP

IIOR = 0
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = ICUNP(ICASE,II1+I)
  IIOR = IIOR+1
  IOCR(IIOR) = JOJ
end do

IINT = IRC(4)
do IQ=1,IINT
  TPQ(IQ) = One
  if (INCPF == 1) TPQ(IQ) = Two/N
  if ((IQ == IREF0) .or. (IP == IREF0)) TPQ(IQ) = Zero
end do
if ((ISDCI == 1) .or. (INCPF == 1) .or. (IP == IREF0)) return

II = 0
IJ = 0
do I=1,LN
  JJ = (IP-1)*LN+I
  if ((ICUNP(ICASE,JJ) == IOCR(I)) .or. (ICUNP(ICASE,JJ) == 3)) cycle
  if (LWSP .and. (ICUNP(ICASE,JJ)*IOCR(I) == 2)) cycle
  if (II == 0) II = I
  IJ = I
end do
NI = IOCR(II)
if (NI > 1) NI = NI-1
NJ = IOCR(IJ)
if (NJ > 1) NJ = NJ-1
do IQ=1,IINT
  IK = 0
  IL = 0
  do I=1,LN
    JJ = (IQ-1)*LN+I
    if ((ICUNP(ICASE,JJ) == IOCR(I)) .or. (ICUNP(ICASE,JJ) == 3)) cycle
    if (LWSP .and. (ICUNP(ICASE,JJ)*IOCR(I) == 2)) cycle
    if (IK == 0) IK = I
    IL = I
  end do
  DIK = Zero
  DIL = Zero
  DJK = Zero
  DJL = Zero
  if (II == IK) DIK = One
  if (II == IL) DIL = One
  if (IJ == IK) DJK = One
  if (IJ == IL) DJL = One
  TPQ(IQ) = (DIK+DIL)/(Two*NI)+(DJK+DJL)/(Two*NJ)
  if (IQ == IREF0) TPQ(IQ) = Zero
end do
if (IPRINT > 15) write(u6,11) (TPQ(IQ),IQ=1,IINT)

return

11 format(5X,'TPQ',10F5.2)

end subroutine TPQSET
