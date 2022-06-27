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

subroutine LOOP16(KM,ISTOP,IT1,IT2)

use guga_global, only: BL1, BL2, BS3, BS4, COUP, COUP1, IA, IB, ICOUP, ICOUP1, IWAY, IY, J1, J2, JM, JM1, K0, K1, K1F, K2, K2F, &
                       K3, K3F
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IDIF, IWAYKM, KM1
real(kind=wp) :: WM0, WP0

ISTOP = 0
KM1 = KM+1
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < 0) .or. (IDIF > 2)) then
  ISTOP = 1
else if (IDIF-1 < 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! NG
    IWAY(KM) = 2
    if ((K1(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)
      COUP(KM) = -COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! PH
    IWAY(KM) = 3
    if ((K3(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K3F(JM(KM1)) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = -BS4(IB(J2(KM1))+1)*COUP(KM1)
    end if
  end if
  if (IWAYKM == 3) ISTOP = 1
else if (IDIF-1 == 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! (OE,QG)
    IWAY(KM) = 2
    if ((K2(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)
      WM0 = BL2(IB(J2(KM1)))
      if (K1F(JM(KM1)) == 0) then
        WP0 = Zero
      else
        WP0 = One/(IB(J2(KM1))+1)
      end if
      COUP(KM) = WM0*COUP1(KM1)+WP0*COUP(KM1)
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)
      WP0 = One/(IB(J2(KM1))+1)
      COUP(KM) = WP0*COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! (RE,JG)
    IWAY(KM) = 3
    if ((K1(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)
      WM0 = -One/(IB(J2(KM1))+1)
      if (K1F(JM(KM1)) == 0) then
        WP0 = Zero
      else
        WP0 = BL1(IB(J2(KM1))+2)
      end if
      COUP(KM) = WM0*COUP1(KM1)+WP0*COUP(KM1)
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)
      WP0 = BL1(IB(J2(KM1))+2)
      COUP(KM) = WP0*COUP(KM1)
    end if
  end if
  if (IWAYKM == 3) then
    ! LH
    IWAY(KM) = 4
    if ((K3(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K3F(JM(KM1)) == 0) then
      IWAYKM = 4
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = -BS4(IB(J2(KM1))+1)*COUP(KM1)
    end if
  end if
  if (IWAYKM == 4) then
    ! PF
    IWAY(KM) = 5
    if ((K3(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else if (K3F(JM1(KM1)) == 0) then
      IWAYKM = 5
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = -BS3(IB(J2(KM1))+1)*COUP1(KM1)
    end if
  end if
  if (IWAYKM == 5) ISTOP = 1
else
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! LF
    IWAY(KM) = 2
    if ((K3(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K3F(JM1(KM1)) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = -BS3(IB(J2(KM1))+1)*COUP1(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! KE
    IWAY(KM) = 3
    if ((K2(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K2F(JM1(KM1)) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)
      COUP(KM) = -COUP1(KM1)
    end if
  end if
  if (IWAYKM == 3) ISTOP = 1
end if

return

end subroutine LOOP16
