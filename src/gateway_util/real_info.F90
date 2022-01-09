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

private
public :: AccMch, ThrInt, PotNuc, Rtrnc, CutInt, TMass, qNuc, PkAcc, &
          Thrs, RadMax, cdMax, EtMax, E1, E2, RPQMin, SadStep, Shake, &
          ChiI2, CoM, rMI, Prin, PAX, CoC, kVector, &
          Real_Info_Dmp, Real_Info_Get

#include "stdalloc.fh"
integer i
real*8 :: AccMch = 1.d-15
real*8 :: ThrInt = 1.d-14
real*8 :: PotNuc = 0.0d0
real*8 :: Rtrnc = 3.d0
real*8 :: CutInt = 1.d-16
real*8 :: TMass = 0.0d0
real*8 :: qNuc = 0.0d0
real*8 :: PkAcc = 1.d-14
real*8 :: Thrs = 1.d-6
real*8 :: RadMax = 0.0d0
real*8 :: cdMax = 0.0d0
real*8 :: EtMax = 0.0d0
real*8 :: E1 = 0.0d0
real*8 :: E2 = 0.0d0
real*8 :: RPQMin = 0.4d0
real*8 :: SadStep = 0.1d0
real*8 :: Shake = -1.0d0
real*8 :: ChiI2 = 0.0d0
real*8 :: CoM(3) = [(0.0d0,i=1,3)]
real*8 :: rMI(6) = [(0.0d0,i=1,6)]
real*8 :: Prin(3) = [(0.0d0,i=1,3)]
real*8 :: PAX(9) = [(0.0d0,i=1,9)]
real*8 :: CoC(3) = [(0.0d0,i=1,3)]
real*8 :: kVector(3) = [(0.0d0,i=1,3)]

contains

subroutine Real_Info_Dmp()

  real*8, allocatable :: rDmp(:)
  integer :: Len = 45

  call mma_allocate(rDmp,Len,Label='rDmp:Real')

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

  call Put_dArray('Real_Info',rDmp,Len)
  call mma_deallocate(rDmp)

end subroutine Real_Info_Dmp

subroutine Real_Info_Get()

  real*8, allocatable :: rDmp(:)
  integer :: Len = 45

  call mma_allocate(rDmp,Len,Label='rDmp:Real')
  call Get_dArray('Real_Info',rDmp,Len)

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
