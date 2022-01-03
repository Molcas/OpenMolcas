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

subroutine LOOP12(KM,ISTOP,IT1,IT2)

use guga_global, only: BS1, BS2, BS3, BS4, COUP, IA, IB, ICOUP, ICOUP1, IPO, IWAY, IY, J1, J2, K0F, K1, K1F, K2, K2F, K3
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IDIF, IWAYKM, J2F, KM1

ISTOP = 0
KM1 = KM+1
J2F = IPO(J2(KM1))
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) then
  ISTOP = 1
else if (IDIF /= 1) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! GC
    IWAY(KM) = 2
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K0F(J2F) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! HD
    IWAY(KM) = 3
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K2F(J2F) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM) = BS4(IB(J2(KM1))+1)*BS2(IB(J2(KM1))+1)*COUP(KM1)
    end if
  end if
  if (IWAYKM == 3) then
    ! GA
    IWAY(KM) = 4
    if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K0F(J2F) == 0) then
      IWAYKM = 4
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = COUP(KM1)
    end if
  end if
  if (IWAYKM == 4) ISTOP = 1
else
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! EC
    IWAY(KM) = 2
    if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K0F(J2F) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! EA
    IWAY(KM) = 3
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K0F(J2F) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = COUP(KM1)
    end if
  end if
  if (IWAYKM == 3) then
    ! FB
    IWAY(KM) = 4
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K1F(J2F) == 0) then
      IWAYKM = 4
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM) = BS3(IB(J2(KM1))+1)*BS1(IB(J2(KM1))+1)*COUP(KM1)
    end if
  end if
  if (IWAYKM == 4) ISTOP = 1
end if

return

end subroutine LOOP12
