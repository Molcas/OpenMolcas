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

subroutine LOOP18(KM,ISTOP,IT1,IT2)

use guga_global, only: BS3, BS4, COUP, COUP1, IA, IB, ICOUP, ICOUP1, IWAY, IY, J1, J2, JM, JM1, K0, K1, K1F, K2, K2F, K3F
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
if ((IDIF < -1) .or. (IDIF > 1)) then
  ISTOP = 1
else if (IDIF < 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! F+H
    IWAY(KM) = 2
    if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K3F(JM(KM1)) == 0) then
      IWAYKM = 2
    else
      J2(KM) = K2(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = BS3(IB(J2(KM1))+3)*BS4(IB(J2(KM1))+1)*COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) ISTOP = 1
else if (IDIF == 0) then
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! (E+E,G+G)
    IWAY(KM) = 2
    if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K2F(JM1(KM1)) /= 0) then
      J2(KM) = K0(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      WM0 = One
      if (K1F(JM(KM1)) == 0) then
        WP0 = Zero
      else
        WP0 = One
      end if
      COUP(KM) = WM0*COUP1(KM1)+WP0*COUP(KM1)
    else if (K1F(JM(KM1)) == 0) then
      IWAYKM = 2
    else
      J2(KM) = K0(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)
      ICOUP(KM) = ICOUP(KM1)
      WP0 = One
      COUP(KM) = WP0*COUP(KM1)
    end if
  end if
  if (IWAYKM == 2) then
    ! F+F
    IWAY(KM) = 3
    if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 3
    else if (K3F(JM1(KM1)) == 0) then
      IWAYKM = 3
    else
      J2(KM) = K1(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = COUP1(KM1)*BS3(IB(J2(KM1))+1)**2
    end if
  end if
  if (IWAYKM == 3) then
    ! H+H
    IWAY(KM) = 4
    if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) then
      IWAYKM = 4
    else if (K3F(JM(KM1)) == 0) then
      IWAYKM = 4
    else
      J2(KM) = K2(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM) = COUP(KM1)*BS4(IB(J2(KM1))+1)**2
    end if
  end if
  if (IWAYKM == 4) ISTOP = 1
else
  IWAYKM = IWAY(KM)
  if (IWAYKM == 1) then
    ! H+F
    IWAY(KM) = 2
    if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) then
      IWAYKM = 2
    else if (K3F(JM1(KM1)) == 0) then
      IWAYKM = 2
    else
      J2(KM) = K1(IT2+J2(KM1))
      J1(KM) = J2(KM)
      ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM) = BS4(IB(J2(KM1))-1)*BS3(IB(J2(KM1))+1)*COUP1(KM1)
    end if
  end if
  if (IWAYKM == 2) ISTOP = 1
end if

return

end subroutine LOOP18
