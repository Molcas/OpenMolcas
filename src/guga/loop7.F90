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

subroutine LOOP7(KM,ISTOP,IT1,IT2)

use guga_global, only: BS1, BS2, COUP, COUP1, IB, ICOUP, ICOUP1, IPO, IWAY, IY, J1, J2, JM, JM1, K0, K1F, K2F, K3
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IWAYKM, J2F, KM1

ISTOP = 0
KM1 = KM+1
J2F = IPO(J2(KM1))
IWAYKM = IWAY(KM)
if (IWAYKM == 1) then
  ! (CB,AD)
  IWAY(KM) = 2
  if ((K0(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) then
    IWAYKM = 2
  else if (K1F(J2F) /= 0) then
    J1(KM) = K0(IT1+J1(KM1))
    J2(KM) = K3(IT2+J2(KM1))
    ICOUP1(KM) = 0
    ICOUP(KM) = IY(IT2+J2(KM1),3)
    COUP1(KM) = BS1(IB(J2(KM1))+1)
    JM1(KM) = K1F(J2F)
    if (K2F(J2F) /= 0) then
      COUP(KM) = BS2(IB(J2(KM1))+1)
      JM(KM) = K2F(J2F)
    end if
  else if (K2F(J2F) == 0) then
    IWAYKM = 2
  else
    J1(KM) = K0(IT1+J1(KM1))
    J2(KM) = K3(IT2+J2(KM1))
    ICOUP1(KM) = 0
    ICOUP(KM) = IY(IT2+J2(KM1),3)
    COUP(KM) = BS2(IB(J2(KM1))+1)
    JM(KM) = K2F(J2F)
  end if
end if
if (IWAYKM == 2) ISTOP = 1

return

end subroutine LOOP7
