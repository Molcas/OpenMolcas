!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module Real_Info

use Constants, only: Zero, One, Three
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: rLen = 45 ! number of elements
real(kind=wp) :: AccMch = 1.0e-15_wp, &
                 cdMax = Zero, &
                 ChiI2 = Zero, &
                 CoC(3) = Zero, &
                 CoM(3) = Zero, &
                 CutInt = 1.0e-16_wp, &
                 E1 = Zero, &
                 E2 = Zero, &
                 EtMax = Zero, &
                 kVector(3) = Zero, &
                 PAX(9) = Zero, &
                 PkAcc = 1.0e-14_wp, &
                 PotNuc = Zero, &
                 Prin(3) = Zero, &
                 qNuc = Zero, &
                 RadMax = Zero, &
                 rMI(6) = Zero, &
                 RPQMin = 0.4_wp, &
                 Rtrnc = Three, &
                 SadStep = 0.1_wp, &
                 Shake = -One, &
                 ThrInt = 1.0e-14_wp, &
                 Thrs = 1.0e-6_wp, &
                 TMass = Zero

public :: AccMch, cdMax, ChiI2, CoC, CoM, CutInt, E1, E2, EtMax, kVector, PAX, PkAcc, PotNuc, Prin, qNuc, RadMax, Real_Info_Dmp, &
          Real_Info_Get, rMI, RPQMin, Rtrnc, SadStep, Shake, ThrInt, Thrs, TMass

contains

subroutine Real_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), allocatable :: rDmp(:)

  call mma_allocate(rDmp,rLen,Label='rDmp:Real')

  rDmp(01) = AccMch
  rDmp(02) = ThrInt
  rDmp(03) = PotNuc
  rDmp(04) = Rtrnc
  rDmp(05) = CutInt
  rDmp(06) = TMass
  rDmp(07) = qNuc
  rDmp(08) = PkAcc
  rDmp(09) = Thrs
  rDmp(10) = RadMax
  rDmp(11) = cdMax
  rDmp(12) = EtMax
  rDmp(13) = E1
  rDmp(14) = E2
  rDmp(15) = RPQMin
  rDmp(16) = SadStep
  rDmp(17) = Shake
  rDmp(18) = ChiI2
  rDmp(19:21) = CoM(1:3)
  rDmp(22:27) = rMI(1:6)
  rDmp(28:30) = Prin(1:3)
  rDmp(31:39) = PAX(1:9)
  rDmp(40:42) = CoC(1:3)
  rDmp(43:45) = kVector(1:3)

  call Put_dArray('Real_Info',rDmp,rLen)
  call mma_deallocate(rDmp)

end subroutine Real_Info_Dmp

subroutine Real_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), allocatable :: rDmp(:)

  call mma_allocate(rDmp,rLen,Label='rDmp:Real')
  call Get_dArray('Real_Info',rDmp,rLen)

  AccMch = rDmp(01)
  ThrInt = rDmp(02)
  PotNuc = rDmp(03)
  Rtrnc = rDmp(04)
  CutInt = rDmp(05)
  TMass = rDmp(06)
  qNuc = rDmp(07)
  PkAcc = rDmp(08)
  Thrs = rDmp(09)
  RadMax = rDmp(10)
  cdMax = rDmp(11)
  EtMax = rDmp(12)
  E1 = rDmp(13)
  E2 = rDmp(14)
  RPQMin = rDmp(15)
  SadStep = rDmp(16)
  Shake = rDmp(17)
  ChiI2 = rDmp(18)
  CoM(1:3) = rDmp(19:21)
  rMI(1:6) = rDmp(22:27)
  Prin(1:3) = rDmp(28:30)
  PAX(1:9) = rDmp(31:39)
  CoC(1:3) = rDmp(40:42)
  kVector(1:3) = rDmp(43:45)

  call mma_deallocate(rDmp)

end subroutine Real_Info_Get

end module Real_Info
