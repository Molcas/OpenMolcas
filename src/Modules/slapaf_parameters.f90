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
Module Slapaf_Parameters
Implicit none
Private
Public:: iRow, iRow_c, iInt, nFix, ddV_Schlegel, HWRS, iOptH, HUpMet, HrmFrq_Show, IRC, &
         nBVec, nDimBC, Curvilinear, Redundant, FindTS, User_Def, Analytic_Hessian, MaxItr, &
         UpMeth, iOptC, HSet, BSet, rHidden, CnstWght, PrQ, lOld, Numerical, Beta, Beta_Disp, &
         Line_Search, iNeg, TSConstraints, GNrm_Threshold, Mode, GrdLbl, GrdMax, &
         StpLbl, StpMax, E_Delta, ThrEne, ThrGrd, nLambda, iRef, ThrCons, ThrMEP, Baker,  &
         eMEPTest, rMEP, MEP, nMEP, MEPNum, MEPCons, dMEPStep, MEP_Type, MEP_Algo, Header, &
         Max_Center, mTROld, Delta, RtRnc, rFuzz, lNmHss, Cubic, lRP, Request_Alaska, Request_RASSI, &
         lOld_Implicit, CallLast, lSoft, lCtoF, Track, TwoRunFiles, isFalcon, Stop

Integer i

Integer:: iRow=0
Integer:: iRow_c=0
Integer:: iInt=0
Integer:: nFix=0
Integer:: IRC=0
Integer:: nBVec=0
Integer:: nDimBC=0
Integer, Parameter:: MaxItr=2000
Integer:: iNeg(2)=[0,0]
Integer:: Mode=-1
Integer:: nLambda=0
Integer:: iRef=0
Integer:: nMEP=MaxItr
Integer:: MEPNum=0
Integer:: Max_Center=15
Integer:: mTROld=0

Logical:: Curvilinear=.True.
Logical:: Redundant=.False.
Logical:: FindTS=.False.
Logical:: HrmFrq_Show=.False.
Logical:: ddV_Schlegel=.False.
Logical:: HWRS=.True.
Logical:: User_Def=.False.
Logical:: Analytic_Hessian=.False.
Logical:: HSet=.False.
Logical:: BSet=.False.
Logical:: PrQ=.False.
Logical:: lOld=.False.
Logical:: Numerical=.False.
Logical:: Line_Search=.True.
Logical:: TSConstraints=.False.
Logical:: Baker=.False.            ! convergence a la Baker
Logical:: eMEPTest=.True.
Logical:: rMEP=.False.
Logical:: MEP=.False.
Logical:: MEPCons=.False.
Logical:: lNmHss=.False.
Logical:: Cubic=.False.
Logical:: lRP=.False.
Logical:: Request_Alaska=.False.
Logical:: Request_RASSI=.False.
Logical:: lOld_Implicit=.False.
Logical:: CallLast=.True.
Logical:: lSoft=.False.
Logical:: lCtoF=.False.
Logical:: Track=.False.
Logical:: TwoRunFiles=.False.
Logical:: isFalcon=.False.
Logical:: Stop=.False.


#include "real.fh"
Real*8:: rHidden=Zero
Real*8:: CnstWght=One
Real*8:: Beta = 0.30D0    !     The threshold for restricted step optimization.
Real*8:: Beta_Disp=0.30D0 !     The threshold for restricted variance optimization.
Real*8:: GNrm_Threshold=0.2D0
Real*8:: GrdMax=Zero, StpMax=Zero
Real*8:: E_Delta=Zero
Real*8:: ThrEne=Zero, ThrGrd=Zero
Real*8:: ThrCons=Zero, ThrMEP=Zero
Real*8:: dMEPStep=0.1D0
Real*8:: Delta=1.0D-2
Real*8:: RtRnc=Three
Real*8:: rFuzz=Half

Character(LEN=8):: GrdLbl='', StpLbl=''
Character(LEN=10):: MEP_TYPE='SPHERE'
Character(LEN=2):: MEP_Algo='GS'
Character(LEN=1):: Header(144)=[('',i=1,144)]
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
!
Integer:: iOptH=4
Character(LEN=6):: HUpMet=' None '
!                                                                      *
!***********************************************************************
!                                                                      *
!.... Optimization method. DO NOT EVER GO BEYOND BIT 30!!!
!
!      iOptC=000000000 (  0) No optimization
!   0  iOptC=000000001 (  1) Quasi Newton-Raphson
!   1  iOptC=000000010 (  2) c1-DIIS
!   2  iOptC=000000100 (  4) c2-DIIS
!   3  iOptC=000001000 (  8) RS-RFO
!   4  iOptC=00001.... ( 16) DIIS, <dx|dx>
!   5  iOptC=00010.... ( 32) DIIS, <dx|g>
!   6  iOptC=00100.... ( 64) DIIS, <g|g>
!   7  iOptC=01....... (128) Minimum, if not set TS search
!   8  iOptC=10....... (256) Optimization with constraint
!   9  iOptC           (512) set: RS-I-RFO, unset: RS-P-RFO
!  10  iOptC          (1024) HMF augmented with weak interactions
!  11  iOptC          (2048) augmented HMF used for selection of
!                            internal coordinates
!  12  iOptC          (4096) set if FindTS
!  13  iOptC          (8192) set if FindTS and in TS regime
!
Integer:: iOptC=2**3 + 2**6 + 2**7 + 2**9 + 2**10 + 2**11
Character(LEN=6):: UpMeth='  RF  '
End Module Slapaf_Parameters
