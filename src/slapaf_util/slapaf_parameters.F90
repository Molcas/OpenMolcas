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

use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp

implicit none
private

! Baker      Convergence a la Baker
! Beta       The threshold for restricted step optimization.
! Beta_Disp  The threshold for restricted variance optimization.
!
!     Hessian update
! 0   iOptH=00000001 (  1) Meyer (disabled)
! 1   iOptH=00000010 (  2) BP (disabled)
! 2   iOptH=00000100 (  4) BFGS
! 3   iOptH=00001000 (  8) None
! 4   iOptH=00010000 ( 16) MPS, for TS search
! 5   iOptH=..1..... ( 32) Not used
! 6   iOptH=01000000 ( 64) EU, for TS search
! 7   iOptH=10000000 (128) TS-BFGS, for TS search
!
! Optimization method. DO NOT EVER GO BEYOND BIT 30!!!
!
!     iOptC=00000000000000  (  0) No optimization
!  0  iOptC=00000000000001  (  1) Quasi Newton-Raphson
!  1  iOptC=00000000000010  (  2) c1-DIIS
!  2  iOptC=00000000000100  (  4) c2-DIIS
!  3  iOptC=00000000001000  (  8) RS-RFO
!  4  iOptC=0000000001....  ( 16) DIIS, <dx|dx>
!  5  iOptC=0000000010....  ( 32) DIIS, <dx|g>
!  6  iOptC=0000000100....  ( 64) DIIS, <g|g>
!  7  iOptC=0000001.......  (128) Minimum, if not set TS search
!  8  iOptC=0000010.......  (256) Optimization with constraint
!  9  iOptC=00001.........  (512) set: RS-I-RFO, unset: RS-P-RFO
! 10  iOptC=0001.......... (1024) HMF augmented with weak interactions
! 11  iOptC=0010.......... (2048) augmented HMF used for selection of internal coordinates
! 12  iOptC=01............ (4096) set if FindTS
! 13  iOptC=10............ (8192) set if FindTS and in TS regime

integer(kind=iwp), parameter :: MaxItr = 2000
integer(kind=iwp) :: iInt = 0, iNeg(2) = 0, iOptC = int(b'111011001000'), iOptH = int(b'100'), IRC = 0, iRef = 0, iRow = 0, &
                     iRow_c = 0, iState(2) = 0, iter = 0, Max_Center = 15, mB_Tot = 0, mdB_Tot = 0, MEPNum = 0, Mode = -1, &
                     mq = 0, mTROld = 0, mTtAtm = 0, MxItr = 0, nBVec = 0, nDimBC = 0, nFix = 0, nLambda = 0, nMEP = MaxItr, &
                     NmIter = 0, nWndw = 5
real(kind=wp) :: Beta = 0.3_wp, Beta_Disp = 0.3_wp, CnstWght = One, Delta = 1.0e-2_wp, dMEPStep = 0.1_wp, E_Delta = Zero, &
                 GNrm_Threshold = 0.2_wp, GrdMax = Zero, StpMax = Zero, rFuzz = Half, rHidden = Zero, RtRnc = Three, &
                 ThrCons = Zero, ThrMEP = Zero, ThrEne = Zero, ThrGrd = Zero
logical(kind=iwp) :: Analytic_Hessian = .false., ApproxNADC = .false., Baker = .false., BSet = .false., CallLast = .true., &
                     Cubic = .false., Curvilinear = .true., ddV_Schlegel = .false., EDiffZero = .false., eMEPTest = .true., &
                     Fallback = .true., FindTS = .false., Force_dB = .false., HrmFrq_Show = .false., HSet = .false., &
                     HWRS = .true., isFalcon = .false., lCtoF = .false., Line_Search = .true., lNmHss = .false., lOld = .false., &
                     lOld_Implicit = .false., lRP = .false., lSoft = .false., MEP = .false., MEPCons = .false., NADC = .false., &
                     Numerical = .false., PrQ = .false., Redundant = .false., Request_Alaska = .false., Request_RASSI = .false., &
                     rMEP = .false., SlStop = .false., Track = .false., TSConstraints = .false., TwoRunFiles = .false., &
                     User_Def = .false., WeightedConstraints = .false.
character(len=10) :: MEP_TYPE = 'SPHERE'
character(len=8) :: GrdLbl = '', StpLbl = ''
character(len=6) :: HUpMet = ' None ', UpMeth = '  RF  '
character(len=2) :: MEP_Algo = 'GS'
character :: Header(144) = ''

public :: Analytic_Hessian, ApproxNADC, Baker, Beta, Beta_Disp, BSet, CallLast, CnstWght, Cubic, Curvilinear, ddV_Schlegel, Delta, &
          dMEPStep, E_Delta, EDiffZero, eMEPTest, Fallback, FindTS, Force_dB, GNrm_Threshold, GrdLbl, GrdMax, Header, HrmFrq_Show, &
          HSet, HUpMet, HWRS, iInt, iNeg, iOptC, iOptH, IRC, iRef, iRow, iRow_c, isFalcon, iState, iter, lCtoF, Line_Search, &
          lNmHss, lOld, lOld_Implicit, lRP, lSoft, Max_Center, MaxItr, mB_Tot, mdB_Tot, MEP, MEP_Algo, MEP_TYPE, MEPCons, MEPNum, &
          Mode, mq, mTROld, mTtAtm, MxItr, NADC, nBVec, nDimBC, nFix, nLambda, nMEP, NmIter, Numerical, nWndw, PrQ, Redundant, &
          Request_Alaska, Request_RASSI, rFuzz, rHidden, rMEP, RtRnc, SlStop, StpLbl, StpMax, ThrCons, ThrEne, ThrGrd, ThrMEP, &
          Track, TSConstraints, TwoRunFiles, UpMeth, User_Def, WeightedConstraints

end module Slapaf_Parameters
