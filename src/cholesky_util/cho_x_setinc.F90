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
                    N1_VecRd, N2_Qual, N2_VecRd, n_MySP, N_Subtr, nBas, nBasT, nChkQ, nCol_BkmThr, nCol_BkmVec, nCol_Chk, nColAB, &
                    nDecom, nDGM_call, nDim_Batch, nInteg, nMisc, nnBstR, nnBstRT, nnShl, nnShl_Tot, nnZTot, nQual_L, nRow_BkmThr, &
                    nRow_BkmVec, nSection, nShell, nSym, nSys_call, NumCho, NumCho_Bak, NumChT, nVec_in_Buf, nVecRS1, RstCho, &
                    RstDia, Run_Mode, ScDiag, ShA, ShAB, ShB, ShC, ShCD, ShD, Span, SSNorm, SSTau, SubScrStat, tDecDrv, tDecom, &
                    Thr_PreScreen, Thr_SimRI, ThrCom, ThrDiag, ThrNeg, TimSec, tInteg, tMisc, Tol_DiaChk, TooNeg, Trace_Idle, &
                    WarNeg, XCho_AdrVec, XDamp, XlDiag, XnBas, XnnShl, XnPass, XnShell, XnSym, XScDiag, XSpan, XThrCom, XThrDiag, &
                    XThrNeg, XTooNeg, XWarNeg
use Constants, only: Zero
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp) :: irc
integer(kind=iwp), parameter :: iLarge = 99999999
real(kind=wp), parameter :: Large = 1.0e15_wp, Small = 1.0e-15_wp

! Set return code.
! ----------------

irc = 0

iPrint = -iLarge

call iZero(iBas,8)
call iZero(nBas,8)
call iZero(XnBas,8)
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

call iZero(LuCho,8)
call iZero(LuSel,8)
call iZero(LuTmp,8)
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
call iZero(iiBstR,8*3)
call iZero(nnBstR,8*3)
call iZero(nnBstRT,3)
mmBstRT = 0
nQual_L(:) = 0
call iZero(iOffQ,8)

call FZero(DiaMax,8)
call FZero(DiaMaxT,8)
call FZero(DiaMin,8)
call FZero(Damp,2)
Span = Large
XlDiag = Large
DiaMnZ = Large
Thr_PreScreen = -Large
iABMnZ = -iLarge
nnZTot = 0

call iZero(NumCho,8)
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
call iZero(iOff_Col,8)

XThrCom = Large
XThrDiag = Large
call FZero(XDamp,2)
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

call iZero(iChkQ,4*(nChkQ+1))
nCol_Chk = -iLarge
call FZero(TimSec,4*nSection)
call FZero(tInteg,2*nInteg)
call FZero(tDecom,2*nDecom)
call FZero(tMisc,2*nMisc)
call FZero(tDecDrv,2)

call iZero(nVecRS1,8)

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

call iZero(ip_ChVBuf_Sym,8)
call iZero(l_ChVBuf_Sym,8)
call iZero(ip_ChVBfI_Sym,8)
call iZero(l_ChVBfI_Sym,8)
call iZero(nVec_in_Buf,8)

Cho_SScreen = .false.
SSTau = Zero
SubScrStat(1) = Zero
SubScrStat(2) = Zero
SSNorm = 'tbp'

call iZero(NumCho_Bak,8)

Cho_Real_Par = .false.

nRow_BkmVec = 0
nCol_BkmVec = 0
nRow_BkmThr = 0
nCol_BkmThr = 0

end subroutine Cho_X_SetInc
