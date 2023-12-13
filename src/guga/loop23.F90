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

subroutine LOOP23(KM,ISTOP,IFAI,IT1,IT2)

use guga_global, only: BL1, BL2, COUP, COUP1, IA, IB, ICOUP, ICOUP1, IWAY, IY, J1, J2, JM, JM1, K0, K0F, K1, K1F, K2, K2F, K3, K3F
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IFAI, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IDIF, IWAYKM, KM1
real(kind=wp) :: WMM, WMP, WPM, WPP

ISTOP = 0
KM1 = KM+1
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < -1) .or. (IDIF > 1)) then
  ISTOP = 1
else if (IDIF < 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! (NQ+,RK+)
    IWAY(KM) = 2
    if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K1F(JM(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM) = -COUP(KM1)/(IB(J2(KM1))+1)
      JM1(KM) = K1F(JM(KM1))
      if (K2F(JM(KM1)) /= 0) then
        COUP(KM) = COUP(KM1)/(IB(J2(KM1))+3)
        JM(KM) = K2F(JM(KM1))
      end if
    else if (K2F(JM(KM1)) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
    end if
  end if
  if (IWAYKM == 2) then
    ! MI+
    IWAY(KM) = 3
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K0F(JM(KM1)) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      COUP(KM) = COUP(KM1)
      JM(KM) = K0F(JM(KM1))
    end if
  end if
  if (IWAYKM == 3) then
    ! NJ+
    IWAY(KM) = 4
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 4
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = -BL1(IB(J2(KM1))+2)*COUP(KM1)
      JM(KM) = K1F(JM(KM1))
    end if
  end if
  if (IWAYKM == 4) then
    ! OK+
    IWAY(KM) = 5
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else if (K2F(JM(KM1)) == 0) then
      IWAYKM = 5
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = -BL2(IB(J2(KM1))+2)*COUP(KM1)
      JM(KM) = K2F(JM(KM1))
    end if
  end if
  if (IWAYKM == 5) then
    ! PL+
    IWAY(KM) = 6
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 6
    else if (K3F(JM(KM1)) == 0) then
      IWAYKM = 6
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM) = COUP(KM1)
      JM(KM) = K3F(JM(KM1))
    end if
  end if
  if (IWAYKM == 6) ISTOP = 1
else if (IDIF == 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! (NN+,RR+JJ+)
    IWAY(KM) = 2
    WMP = Zero
    WPP = Zero
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K1F(JM1(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM) = COUP1(KM1)
      JM1(KM) = K1F(JM1(KM1))
      if (K2F(JM1(KM1)) /= 0) then
        WMP = One/((IB(J2(KM1))+1)**2)
        JM(KM) = K2F(JM1(KM1))
        if (K1F(JM(KM1)) /= 0) then
          WPP = BL1(IB(J2(KM1))+2)**2
          JM(KM) = K1F(JM(KM1))
        end if
        COUP(KM) = WMP*COUP1(KM1)+WPP*COUP(KM1)
      else if (K1F(JM(KM1)) /= 0) then
        WPP = BL1(IB(J2(KM1))+2)**2
        JM(KM) = K1F(JM(KM1))
        COUP(KM) = WPP*COUP(KM1)
      end if
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WMP = One/((IB(J2(KM1))+1)**2)
      JM(KM) = K2F(JM1(KM1))
      if (K1F(JM(KM1)) /= 0) then
        WPP = BL1(IB(J2(KM1))+2)**2
        JM(KM) = K1F(JM(KM1))
      end if
      COUP(KM) = WMP*COUP1(KM1)+WPP*COUP(KM1)
    else if (K1F(JM(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WPP = BL1(IB(J2(KM1))+2)**2
      JM(KM) = K1F(JM(KM1))
      COUP(KM) = WPP*COUP(KM1)
    else if (IFAI == 0) then
      IWAYKM = 2
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM) = Zero
      COUP(KM) = Zero
    end if
  end if
  if (IWAYKM == 2) then
    ! (KK+,OO+,QQ+)
    IWAY(KM) = 3
    WMM = Zero
    WPM = Zero
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K2F(JM(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = COUP(KM1)
      JM(KM) = K2F(JM(KM1))
      if (K2F(JM1(KM1)) /= 0) then
        WMM = BL2(IB(J2(KM1)))**2
        JM1(KM) = K2F(JM1(KM1))
        if (K1F(JM(KM1)) /= 0) then
          WPM = One/((IB(J2(KM1))+1)**2)
          JM1(KM) = K1F(JM(KM1))
        end if
        COUP1(KM) = WMM*COUP1(KM1)+WPM*COUP(KM1)
      else if (K1F(JM(KM1)) /= 0) then
        WPM = One/((IB(J2(KM1))+1)**2)
        JM1(KM) = K1F(JM(KM1))
        COUP1(KM) = WPM*COUP(KM1)
      end if
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WMM = BL2(IB(J2(KM1)))**2
      JM1(KM) = K2F(JM1(KM1))
      if (K1F(JM(KM1)) /= 0) then
        WPM = One/((IB(J2(KM1))+1)**2)
        JM1(KM) = K1F(JM(KM1))
      end if
      COUP1(KM) = WMM*COUP1(KM1)+WPM*COUP(KM1)
    else if (K1F(JM(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WPM = One/((IB(J2(KM1))+1)**2)
      JM1(KM) = K1F(JM(KM1))
      COUP1(KM) = WPM*COUP(KM1)
    else if (IFAI == 0) then
      IWAYKM = 3
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM) = Zero
      COUP(KM) = Zero
    end if
  end if
  if (IWAYKM == 3) then
    ! (MM+,II+)
    IWAY(KM) = 4
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K0F(JM1(KM1)) /= 0) then
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      COUP1(KM) = COUP1(KM1)
      JM1(KM) = K0F(JM1(KM1))
      if (K0F(JM(KM1)) /= 0) then
        COUP(KM) = COUP(KM1)
        JM(KM) = K0F(JM(KM1))
      end if
    else if (K0F(JM(KM1)) /= 0) then
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      COUP(KM) = COUP(KM1)
      JM(KM) = K0F(JM(KM1))
    else if (IFAI == 0) then
      IWAYKM = 4
    else
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      COUP1(KM) = Zero
      COUP(KM) = Zero
    end if
  end if
  if (IWAYKM == 4) then
    ! (PP+,LL+)
    IWAY(KM) = 5
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else if (K3F(JM1(KM1)) /= 0) then
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM) = COUP1(KM1)
      JM1(KM) = K3F(JM1(KM1))
      if (K3F(JM(KM1)) /= 0) then
        COUP(KM) = COUP(KM1)
        JM(KM) = K3F(JM(KM1))
      end if
    else if (K3F(JM(KM1)) /= 0) then
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM) = COUP(KM1)
      JM(KM) = K3F(JM(KM1))
    else if (IFAI == 0) then
      IWAYKM = 5
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM) = Zero
      COUP(KM) = Zero
    end if
  end if
  if (IWAYKM == 5) then
    ! (OR+,QJ+)
    IWAY(KM) = 6
    if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 6
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WMP = -BL2(IB(J2(KM1)))/(IB(J2(KM1))+1)
      if (K1F(JM(KM1)) == 0) then
        WPP = Zero
        JM(KM) = K2F(JM1(KM1))
      else
        WPP = BL1(IB(J2(KM1))+2)/(IB(J2(KM1))+1)
        JM(KM) = K1F(JM(KM1))
      end if
      COUP(KM) = WMP*COUP1(KM1)+WPP*COUP(KM1)
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 6
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WPP = BL1(IB(J2(KM1))+2)/(IB(J2(KM1))+1)
      JM(KM) = K1F(JM(KM1))
      COUP(KM) = WPP*COUP(KM1)
    end if
  end if
  if (IWAYKM == 6) then
    ! (RO+,JQ+)
    IWAY(KM) = 7
    if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 7
    else if (K2F(JM1(KM1)) /= 0) then
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WMM = -BL2(IB(J2(KM1)))/(IB(J2(KM1))+1)
      if (K1F(JM(KM1)) == 0) then
        WPM = Zero
        JM1(KM) = K2F(JM1(KM1))
      else
        WPM = BL1(IB(J2(KM1))+2)/(IB(J2(KM1))+1)
        JM1(KM) = K1F(JM(KM1))
      end if
      COUP1(KM) = WMM*COUP1(KM1)+WPM*COUP(KM1)
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 7
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WPM = BL1(IB(J2(KM1))+2)/(IB(J2(KM1))+1)
      JM1(KM) = K1F(JM(KM1))
      COUP1(KM) = WPM*COUP(KM1)
    end if
  end if
  if (IWAYKM == 7) ISTOP = 1
else
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! (QN+,KR+)
    IWAY(KM) = 2
    if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K1F(JM1(KM1)) /= 0) then
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM) = -COUP1(KM1)/(IB(J2(KM1))-1)
      JM1(KM) = K1F(JM1(KM1))
      if (K2F(JM1(KM1)) /= 0) then
        COUP(KM) = COUP1(KM1)/(IB(J2(KM1))+1)
        JM(KM) = K2F(JM1(KM1))
      end if
    else if (K2F(JM1(KM1)) == 0) then
      IWAYKM = 2
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = COUP1(KM1)/(IB(J2(KM1))+1)
      JM(KM) = K2F(JM1(KM1))
    end if
  end if
  if (IWAYKM == 2) then
    ! IM+
    IWAY(KM) = 3
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K0F(JM1(KM1)) == 0) then
      IWAYKM = 3
    else
      J1(KM) = K0(IT1+J1(KM1))
      J2(KM) = K0(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      COUP1(KM) = COUP1(KM1)
      JM1(KM) = K0F(JM1(KM1))
    end if
  end if
  if (IWAYKM == 3) then
    ! JN+
    IWAY(KM) = 4
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K1F(JM1(KM1)) == 0) then
      IWAYKM = 4
    else
      J1(KM) = K1(IT1+J1(KM1))
      J2(KM) = K1(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM) = -BL1(IB(J2(KM1)))*COUP1(KM1)
      JM1(KM) = K1F(JM1(KM1))
    end if
  end if
  if (IWAYKM == 4) then
    ! KO+
    IWAY(KM) = 5
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 5
    else if (K2F(JM1(KM1)) == 0) then
      IWAYKM = 5
    else
      J1(KM) = K2(IT1+J1(KM1))
      J2(KM) = K2(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM) = -BL2(IB(J2(KM1)))*COUP1(KM1)
      JM1(KM) = K2F(JM1(KM1))
    end if
  end if
  if (IWAYKM == 5) then
    ! LP+
    IWAY(KM) = 6
    if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
      IWAYKM = 6
    else if (K3F(JM1(KM1)) == 0) then
      IWAYKM = 6
    else
      J1(KM) = K3(IT1+J1(KM1))
      J2(KM) = K3(IT2+J2(KM1))
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM) = COUP1(KM1)
      JM1(KM) = K3F(JM1(KM1))
    end if
  end if
  if (IWAYKM == 6) ISTOP = 1
end if

return

end subroutine LOOP23
