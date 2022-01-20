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
!
! Some incomplete documentation (logical variables)
!
! Vlct       : logical flag to indicate that the velocity integrals will
!              be computed
! lRel       : logical flag to indicate that the 1-electron Darwin and
!              mass-velocity integrals should be computed
! UnNorm     : logical flag to indicate that the primitive basis
!              functions should not be normalized
! lSchw      : logical flag to indicate that the Cauchy-Schwartz
!              inequality should be used in the prescreening
! lAMFI      : logical flag for Atomic Mean Field Integrals
! GIAO       : integrals for first derivative with respect to B
! Do_FckInt  : logical flag for FckInt
! Do_GuessOrb: logical flag for GuessOrb
! EMFR       : electromagnetic field radiation
! FNMC       : finite nuclear mass correction

module Gateway_Info

use Constants, only: Zero, One, Three
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: lLen = 19, & ! number of logical elements
                                rLen = 45 ! number of real, elements
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
                 PAX(3,3) = Zero, &
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
logical(kind=iwp) :: Align_Only = .false., &
                     Do_Align = .true., &
                     Do_FckInt = .true., &
                     Do_GuessOrb = .true., &
                     DoFMM = .false., &
                     EMFR = .false., &
                     FNMC = .false., &
                     GIAO = .false., &
                     lAMFI = .false., &
                     lDOWNONLY = .false., &
                     lMXTC = .false., &
                     lRel = .false., &
                     lRP = .false., &
                     lRP_Post = .false., &
                     lSchw = .true., &
                     lUPONLY = .false., &
                     NEMO = .false., &
                     UnNorm = .false., &
                     Vlct = .true.

public :: AccMch, Align_Only, cdMax, ChiI2, CoC, CoM, CutInt, Do_Align, Do_FckInt, Do_GuessOrb, DoFMM, E1, E2, EMFR, EtMax, FNMC, &
          GIAO, kVector, lAMFI, lDOWNONLY, lMXTC, lRel, lRP, lRP_Post, lSchw, lUPONLY, NEMO, PAX, PkAcc, PotNuc, Prin, qNuc, &
          RadMax, Gateway_Info_Dmp, Gateway_Info_Get, rMI, RPQMin, Rtrnc, SadStep, Shake, ThrInt, Thrs, TMass, UnNorm, Vlct

contains

subroutine Gateway_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), allocatable :: iDmp(:)
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
  rDmp(31:39) = reshape(PAX(1:3,1:3),[3*3])
  rDmp(40:42) = CoC(1:3)
  rDmp(43:45) = kVector(1:3)

  call Put_dArray('Real_Info',rDmp,rLen)
  call mma_deallocate(rDmp)

  call mma_allocate(iDmp,lLen,Label='iDmp:Logical')

  iDmp(01) = merge(1,0,Vlct)
  iDmp(02) = merge(1,0,lRel)
  iDmp(03) = merge(1,0,UnNorm)
  iDmp(04) = merge(1,0,lSchw)
  iDmp(05) = merge(1,0,lAMFI)
  iDmp(06) = merge(1,0,NEMO)
  iDmp(07) = merge(1,0,Do_GuessOrb)
  iDmp(08) = merge(1,0,Do_FckInt)
  iDmp(09) = merge(1,0,Align_Only)
  iDmp(10) = merge(1,0,DoFMM)
  iDmp(11) = merge(1,0,lRP)
  iDmp(12) = merge(1,0,lRP_Post)
  iDmp(13) = merge(1,0,Do_Align)
  iDmp(14) = merge(1,0,EMFR)
  iDmp(15) = merge(1,0,GIAO)
  iDmp(16) = merge(1,0,lUPONLY)
  iDmp(17) = merge(1,0,lDOWNONLY)
  iDmp(18) = merge(1,0,FNMC)
  iDmp(19) = merge(1,0,lMXTC)

  call Put_iArray('Logical_Info',iDmp,lLen)
  call mma_deallocate(iDmp)

end subroutine Gateway_Info_Dmp

subroutine Gateway_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), allocatable :: iDmp(:)
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
  PAX(1:3,1:3) = reshape(rDmp(31:39),[3,3])
  CoC(1:3) = rDmp(40:42)
  kVector(1:3) = rDmp(43:45)

  call mma_deallocate(rDmp)

  call mma_allocate(iDmp,lLen,Label='iDmp:Logical')
  call Get_iArray('Logical_Info',iDmp,lLen)

  Vlct = iDmp(1) > 0
  lRel = iDmp(2) > 0
  UnNorm = iDmp(3) > 0
  lSchw = iDmp(4) > 0
  lAMFI = iDmp(5) > 0
  NEMO = iDmp(6) > 0
  Do_GuessOrb = iDmp(7) > 0
  Do_FckInt = iDmp(8) > 0
  Align_Only = iDmp(9) > 0
  DoFMM = iDmp(10) > 0
  lRP = iDmp(11) > 0
  lRP_Post = iDmp(12) > 0
  Do_Align = iDmp(13) > 0
  EMFR = iDmp(14) > 0
  GIAO = iDmp(15) > 0
  lUPONLY = iDmp(16) > 0
  lDOWNONLY = iDmp(17) > 0
  FNMC = iDmp(18) > 0
  lMXTC = iDmp(19) > 0

  call mma_deallocate(iDmp)

end subroutine Gateway_Info_Get

end module Gateway_Info
