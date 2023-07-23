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
! Purpose: define all entries in include files
!          choprint.fh
!          choorb.fh
!          cholesky.fh
!          chosew.fh
!          chovecbuf.f90
!          chosubscr.fh
!          chpari.fh
!          cho_para_info.fh
!          and some in the Module choarr.f90

use ChoArr, only: n_MySP, nDim_Batch, nQual_L
use ChoBkm, only: nCol_BkmThr, nCol_BkmVec, nRow_BkmThr, nRow_BkmVec
use ChoSubScr, only: Cho_SScreen, SSNorm, SSTau, SubScrStat
use ChoVecBuf, only: ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, l_CHVBUF_SYM, nVec_in_Buf
use ChPari, only: NumCho_Bak
use Constants, only: Zero
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp) :: irc
#include "choorb.fh"
#include "choprint.fh"
#include "cholesky.fh"
#include "cho_para_info.fh"
integer(kind=iwp), parameter :: iLarge = 99999999
real(kind=wp), parameter :: Large = 1.0e15_wp, Small = 1.0e-15_wp

! Set return code.
! ----------------

irc = 0

! choprint.fh.
! -------------

iPrint = -iLarge

! choorb.fh.
! -----------

call iZero(iBas,8)
call iZero(nBas,8)
call iZero(XnBas,8)
nBasT = 0

! cholesky.fh.
! -------------

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

! chovecbuf.f90.
! --------------

call iZero(ip_ChVBuf_Sym,8)
call iZero(l_ChVBuf_Sym,8)
call iZero(ip_ChVBfI_Sym,8)
call iZero(l_ChVBfI_Sym,8)
call iZero(nVec_in_Buf,8)

! chosubscr.fh.
! -------------

Cho_SScreen = .false.
SSTau = Zero
SubScrStat(1) = Zero
SubScrStat(2) = Zero
SSNorm = 'tbp'

! Module chpari
! -------------

call iZero(NumCho_Bak,8)

! cho_para_info.fh.
! -----------------

Cho_Real_Par = .false.

! chobkm.f90
! ----------

nRow_BkmVec = 0
nCol_BkmVec = 0
nRow_BkmThr = 0
nCol_BkmThr = 0

end subroutine Cho_X_SetInc
