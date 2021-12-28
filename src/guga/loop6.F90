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

subroutine LOOP6(KM,ISTOP,IT1,IT2)

use guga_global, only: BL1, BL2, COUP, IA, IB, ICOUP, ICOUP1, IWAY, IY, J1, J2, K0, K1, K2, K3
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IDIF, IWAYKM, KM1

ISTOP = 0
KM1 = KM+1
IDIF = IA(J2(KM1))-IA(J1(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) then
  ISTOP = 1
else if (IDIF /= 0) then
  ! CASE I-L AND Q
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    IWAY(KM) = 2
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else
      COUP(KM) = COUP(KM1)
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    IWAY(KM) = 3
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else
      COUP(KM) = BL1(IB(J1(KM1))+1)*COUP(KM1)
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
    end if
  end if
  if (IWAYKM == 3) then
    IWAY(KM) = 4
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else
      COUP(KM) = -COUP(KM1)
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
    end if
  end if
  if (IWAYKM == 4) then
    IWAY(KM) = 5
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else
      COUP(KM) = -COUP(KM1)
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
    end if
  end if
  if (IWAYKM == 5) then
    IWAY(KM) = 6
    if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 6
    else
      COUP(KM) = COUP(KM1)/IB(J1(KM1))
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
    end if
  end if
  if (IWAYKM == 6) ISTOP = 1
else
  ! CASE M-P AND R
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    IWAY(KM) = 2
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else
      COUP(KM) = COUP(KM1)
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    IWAY(KM) = 3
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else
      COUP(KM) = -COUP(KM1)
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
    end if
  end if
  if (IWAYKM == 3) then
    IWAY(KM) = 4
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else
      COUP(KM) = BL2(IB(J1(KM1))+1)*COUP(KM1)
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
    end if
  end if
  if (IWAYKM == 4) then
    IWAY(KM) = 5
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else
      COUP(KM) = -COUP(KM1)
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
    end if
  end if
  if (IWAYKM == 5) then
    IWAY(KM) = 6
    if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 6
    else
      COUP(KM) = -COUP(KM1)/(IB(J1(KM1))+2)
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
    end if
  end if
  if (IWAYKM == 6) ISTOP = 1
end if

return

end subroutine LOOP6
