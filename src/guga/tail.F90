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
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine TAIL(LL,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)

use guga_global, only: ICOUP, ICOUP1, IWAY, IY, J2, LN
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LL, IJJ, ITAIL, L0(*), L1(*), L2(*), L3(*), IT1, IT2
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: I, IWAYKM, KM, KM1, L
logical(kind=iwp) :: first

L = LL+1
if (ITAIL == 0) then
  return
end if
do I=1,ITAIL
  ITAI(I) = 0
end do
if (L == LN+1) ITAI(1) = 1
KM = L
J2(KM) = IJJ
ICOUP(L) = 1
ICOUP1(L) = 1
first = .true.
do
  if (first) then
    KM = KM+1
    IWAY(KM) = 0
    first = .false.
  end if
  KM1 = KM-1
  IWAYKM = IWAY(KM)
  if (IWAYKM < 1) then
    if ((L0(IT1+J2(KM1)) == 0) .or. (L0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 1
    else
      J2(KM) = L0(IT2+J2(KM1))
      IWAY(KM) = 1
      ICOUP(KM) = ICOUP(KM1)
      ICOUP1(KM) = ICOUP1(KM1)
    end if
  end if
  if (IWAYKM == 1) then
    if ((L1(IT1+J2(KM1)) == 0) .or. (L1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else
      J2(KM) = L1(IT2+J2(KM1))
      IWAY(KM) = 2
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),1)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),1)
    end if
  end if
  if (IWAYKM == 2) then
    if ((L2(IT1+J2(KM1)) == 0) .or. (L2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else
      J2(KM) = L2(IT2+J2(KM1))
      IWAY(KM) = 3
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),2)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),2)
    end if
  end if
  if (IWAYKM == 3) then
    if ((L3(IT1+J2(KM1)) == 0) .or. (L3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else
      J2(KM) = L3(IT2+J2(KM1))
      IWAY(KM) = 4
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),3)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),3)
    end if
  end if
  if (IWAYKM > 3) then
    KM = KM-1
    if (KM == L) exit
  else
    if (KM /= LN+1) then
      first = .true.
    else
      ITAI(ICOUP(LN+1)) = ICOUP1(LN+1)
    end if
  end if
end do

return

end subroutine TAIL
