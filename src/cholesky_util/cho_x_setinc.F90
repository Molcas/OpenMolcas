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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_X_SetInc(irc)
!
! T.B. Pedersen, July 2004.
!
! Purpose: define some entries in the module Cholesky

use Cholesky, only: BlockSize, ChkOnly, Cho_1Center, Cho_AdrVec, Cho_DecAlg, Cho_DecAlg_Def, Cho_DiaChk, Cho_Fake_Par, Cho_IntChk, &
                    Cho_IOVec, Cho_MinChk, Cho_No2Center, Cho_PreScreen, Cho_Real_Par, Cho_ReOrd, Cho_SimP, Cho_SimRI, &
                    Cho_SScreen, Cho_TrcNeg, Cho_TstScreen, Cho_UseAbs, Damp, DiaMax, DiaMaxT, DiaMin, DiaMin, DiaMnZ, Did_DecDrv, &
                    Frac_ChVBuf, Haltit, iABMnZ, iAlQua, iBas, iChkQ, IFCSew, iiBstR, iOff_Col, iOffQ, ip_CHVBFI_SYM, &
                    ip_CHVBUF_SYM, IPRINT, l_CHVBFI_SYM, l_CHVBUF_SYM, lBuf, LuCho, LuMap, LuPri, LuRed, LuRst, LuScr, LuSel, &
                    LuTmp, MaxQual, MaxRed, MaxVec, MinQual, mmBstRT, Mode_Screen, ModRst, Mx2Sh, MxORSh, MxShPr, N1_Qual, &
                    N1_VecRd, N2_Qual, N2_VecRd, n_MySP, N_Subtr, nBas, nBasT, nCol_BkmThr, nCol_BkmVec, nCol_Chk, nColAB, &
                    nDGM_call, nDim_Batch, nnBstR, nnBstRT, nnShl, nnShl_Tot, nnZTot, nQual_L, nRow_BkmThr, nRow_BkmVec, nShell, &
                    nSym, nSys_call, NumCho, NumCho_Bak, NumChT, nVec_in_Buf, nVecRS1, RstCho, RstDia, Run_Mode, ScDiag, ShA, &
                    ShAB, ShB, ShC, ShCD, ShD, Span, SSNorm, SSTau, SubScrStat, tDecDrv, tDecom, Thr_PreScreen, Thr_SimRI, ThrCom, &
                    ThrDiag, ThrNeg, TimSec, tInteg, tMisc, Tol_DiaChk, TooNeg, Trace_Idle, WarNeg, XCho_AdrVec, XDamp, XlDiag, &
                    XnBas, XnnShl, XnPass, XnShell, XnSym, XScDiag, XSpan, XThrCom, XThrDiag, XThrNeg, XTooNeg, XWarNeg
use Constants, only: Zero
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), parameter :: iLarge = 99999999
real(kind=wp), parameter :: Large = 1.0e15_wp, Small = 1.0e-15_wp

! Set return code.
! ----------------

irc = 0

iPrint = -iLarge

iBas(:) = 0
nBas(:) = 0
XnBas(:) = 0
nBasT = 0

ThrCom = Large
ThrDiag = Large
Tol_DiaChk = -Large
ThrNeg = Large
WarNeg = Large
TooNeg = Large
nSym = -iLarge
lBuf = -iLarge
MinQual = iLarge
MaxQUal = iLarge
IFCSew = -iLarge
Mode_Screen = -iLarge
Cho_DecAlg = -iLarge
Cho_DecAlg_Def = -iLarge
iAlQua = iLarge
MxShPr = iLarge
ModRst = iLarge
Run_Mode = iLarge
ScDiag = .false.
ChkOnly = .false.
Cho_IntChk = .false.
Cho_MinChk = .false.
Cho_UseAbs = .false.
Cho_TrcNeg = .false.
Cho_ReOrd = .false.
Cho_DiaChk = .false.
Cho_TstScreen = .false.
Cho_1Center = .false.
Cho_No2Center = .false.
Cho_PreScreen = .false.
Cho_SimP = .false.
Cho_Fake_Par = .false.
RstDia = .false.
RstCho = .false.
Did_DecDrv = .false.
HaltIt = .false.
Trace_Idle = .false.

LuCho(:) = 0
LuSel(:) = 0
LuTmp(:) = 0
LuPri = 0
LuScr = 0
LuRed = 0
LuRst = 0
LuMap = 0

nShell = 0
nnShl_Tot = 0
nnShl = 0
MxORSh = 0
Mx2Sh = 0
iiBstR(:,:) = 0
nnBstR(:,:) = 0
nnBstRT(:) = 0
mmBstRT = 0
nQual_L(:) = 0
iOffQ(:) = 0

DiaMax(:) = Zero
DiaMaxT(:) = Zero
DiaMin(:) = Zero
Damp(:) = Zero
Span = Large
XlDiag = Large
DiaMnZ = Large
Thr_PreScreen = -Large
iABMnZ = -iLarge
nnZTot = 0

NumCho(:) = 0
NumChT = 0
MaxVec = 0
MaxRed = 0
BlockSize = -iLarge

ShA = -iLarge
ShB = -iLarge
ShC = -iLarge
ShD = -iLarge
ShAB = -iLarge
ShCD = -iLarge
nColAB = -iLarge
iOff_Col(:) = 0

XThrCom = Large
XThrDiag = Large
XDamp(:) = Zero
XSpan = Large
XThrNeg = Large
XWarNeg = Large
XTooNeg = Large
XnSym = 0
XnShell = 0
XnnShl = 0
XnPass = 0
XScDiag = .false.
XCho_AdrVec = -iLarge

iChkQ(:,:) = 0
nCol_Chk = -iLarge
TimSec(:,:) = Zero
tInteg(:,:) = Zero
tDecom(:,:) = Zero
tMisc(:,:) = Zero
tDecDrv(:) = Zero

nVecRS1(:) = 0

Cho_AdrVec = -iLarge
Cho_IOVec = -iLarge
nSys_Call = 0
nDGM_Call = 0
N1_VecRd = 0
N2_VecRd = 0
N_Subtr = 0

N1_Qual = -iLarge
N2_Qual = iLarge

Frac_ChVBuf = Zero

nDim_Batch(:) = 0

nQual_L(:) = 0

n_MySP = 0

Cho_SimRI = .false.
Thr_SimRI = -Large

ip_ChVBuf_Sym(:) = 0
l_ChVBuf_Sym(:) = 0
ip_ChVBfI_Sym(:) = 0
l_ChVBfI_Sym(:) = 0
nVec_in_Buf(:) = 0

Cho_SScreen = .false.
SSTau = Zero
SubScrStat(1) = Zero
SubScrStat(2) = Zero
SSNorm = 'tbp'

NumCho_Bak(:) = 0

Cho_Real_Par = .false.

nRow_BkmVec = 0
nCol_BkmVec = 0
nRow_BkmThr = 0
nCol_BkmThr = 0

end subroutine Cho_X_SetInc
