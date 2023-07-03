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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module Slapaf_Parameters

implicit none
private

public :: iRow, iRow_c, iInt, nFix, ddV_Schlegel, HWRS, iOptH, HUpMet, HrmFrq_Show, IRC, nBVec, nDimBC, Curvilinear, Redundant, &
          FindTS, User_Def, Analytic_Hessian, MaxItr, UpMeth, iOptC, HSet, BSet, rHidden, CnstWght, PrQ, lOld, Numerical, Beta, &
          Beta_Disp, Line_Search, iNeg, TSConstraints, GNrm_Threshold, Mode, GrdLbl, GrdMax, StpLbl, StpMax, E_Delta, ThrEne, &
          ThrGrd, nLambda, iRef, ThrCons, ThrMEP, Baker, eMEPTest, rMEP, MEP, nMEP, MEPNum, MEPCons, dMEPStep, MEP_Type, MEP_Algo, &
          Header, Max_Center, mTROld, Delta, RtRnc, rFuzz, lNmHss, Cubic, lRP, Request_Alaska, Request_RASSI, lOld_Implicit, &
          CallLast, lSoft, lCtoF, Track, TwoRunFiles, isFalcon, stop, NmIter, MxItr, mTtAtm, nWndw, iter, WeightedConstraints, &
          mB_Tot, mdB_Tot, mq, Force_dB, NADC, EDiffZero, ApproxNADC, iState, Fallback

integer :: i

integer :: iRow = 0
integer :: iRow_c = 0
integer :: iInt = 0
integer :: nFix = 0
integer :: IRC = 0
integer :: nBVec = 0
integer :: nDimBC = 0
integer, parameter :: MaxItr = 2000
integer :: iNeg(2) = [0,0]
integer :: Mode = -1
integer :: nLambda = 0
integer :: iRef = 0
integer :: nMEP = MaxItr
integer :: MEPNum = 0
integer :: Max_Center = 15
integer :: mTROld = 0
integer :: NmIter = 0
integer :: MxItr = 0
integer :: mTtAtm = 0
integer :: nWndw = 5
integer :: iter = 0
integer :: mB_Tot = 0
integer :: mdB_Tot = 0
integer :: mq = 0
integer :: iState(2) = [0,0]

logical :: Curvilinear = .true.
logical :: Redundant = .false.
logical :: FindTS = .false.
logical :: HrmFrq_Show = .false.
logical :: ddV_Schlegel = .false.
logical :: HWRS = .true.
logical :: User_Def = .false.
logical :: Analytic_Hessian = .false.
logical :: HSet = .false.
logical :: BSet = .false.
logical :: PrQ = .false.
logical :: lOld = .false.
logical :: Numerical = .false.
logical :: Line_Search = .true.
logical :: TSConstraints = .false.
logical :: Baker = .false.            ! convergence a la Baker
logical :: eMEPTest = .true.
logical :: rMEP = .false.
logical :: MEP = .false.
logical :: MEPCons = .false.
logical :: lNmHss = .false.
logical :: Cubic = .false.
logical :: lRP = .false.
logical :: Request_Alaska = .false.
logical :: Request_RASSI = .false.
logical :: lOld_Implicit = .false.
logical :: CallLast = .true.
logical :: lSoft = .false.
logical :: lCtoF = .false.
logical :: Track = .false.
logical :: TwoRunFiles = .false.
logical :: isFalcon = .false.
logical :: stop = .false.
logical :: WeightedConstraints = .false.
logical :: Force_dB = .false.
logical :: NADC = .false.
logical :: EDiffZero = .false.
logical :: ApproxNADC = .false.
logical :: Fallback = .true.

#include "real.fh"
real*8 :: rHidden = Zero
real*8 :: CnstWght = One
real*8 :: Beta = 0.30d0      !     The threshold for restricted step optimization.
real*8 :: Beta_Disp = 0.30d0 !     The threshold for restricted variance optimization.
real*8 :: GNrm_Threshold = 0.2d0
real*8 :: GrdMax = Zero, StpMax = Zero
real*8 :: E_Delta = Zero
real*8 :: ThrEne = Zero, ThrGrd = Zero
real*8 :: ThrCons = Zero, ThrMEP = Zero
real*8 :: dMEPStep = 0.1d0
real*8 :: Delta = 1.0D-2
real*8 :: RtRnc = Three
real*8 :: rFuzz = Half

character(len=8) :: GrdLbl = '', StpLbl = ''
character(len=10) :: MEP_TYPE = 'SPHERE'
character(len=2) :: MEP_Algo = 'GS'
character(len=1) :: Header(144) = [('',i=1,144)]
!                                                                      *
!***********************************************************************
!                                                                      *
!     Hessian update
! 1   iOptH=00000001 (  1) Meyer (disabled)
! 2   iOptH=00000010 (  2) BP (disabled)
! 3   iOptH=00000100 (  4) BFGS
! 4   iOptH=00001000 (  8) None
! 5   iOptH=00010000 ( 16) MPS, for TS search
! 6   iOptH=-.1..... ( 32) Not used
! 7   iOptH=01000000 ( 64) EU, for TS search
! 8   iOptH=10000000 (128) TS-BFGS, for TS search

integer :: iOptH = 4
character(len=6) :: HUpMet = ' None '
!                                                                      *
!***********************************************************************
!                                                                      *
! Optimization method. DO NOT EVER GO BEYOND BIT 30!!!
!
!     iOptC=000000000 (  0) No optimization
!  0  iOptC=000000001 (  1) Quasi Newton-Raphson
!  1  iOptC=000000010 (  2) c1-DIIS
!  2  iOptC=000000100 (  4) c2-DIIS
!  3  iOptC=000001000 (  8) RS-RFO
!  4  iOptC=00001.... ( 16) DIIS, <dx|dx>
!  5  iOptC=00010.... ( 32) DIIS, <dx|g>
!  6  iOptC=00100.... ( 64) DIIS, <g|g>
!  7  iOptC=01....... (128) Minimum, if not set TS search
!  8  iOptC=10....... (256) Optimization with constraint
!  9  iOptC           (512) set: RS-I-RFO, unset: RS-P-RFO
! 10  iOptC          (1024) HMF augmented with weak interactions
! 11  iOptC          (2048) augmented HMF used for selection of internal coordinates
! 12  iOptC          (4096) set if FindTS
! 13  iOptC          (8192) set if FindTS and in TS regime

integer :: iOptC = 2**3+2**6+2**7+2**9+2**10+2**11
character(len=6) :: UpMeth = '  RF  '

end module Slapaf_Parameters
