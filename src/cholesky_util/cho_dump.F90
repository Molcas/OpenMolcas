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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_Dump(irc,Lunit)
!
! T.B. Pedersen, March 2005.
!
! Purpose: print some entries in the module Cholesky
!
! On input, Lunit is the logical unit to print to...

use Cholesky, only: BlockSize, ChkOnly, Cho_1Center, Cho_AdrVec, Cho_DecAlg, Cho_DecAlg_Def, Cho_DiaChk, Cho_Fake_Par, Cho_IntChk, &
                    Cho_IOVec, Cho_MinChk, Cho_No2Center, Cho_PreScreen, Cho_Reord, Cho_SimP, Cho_SScreen, Cho_TrcNeg, &
                    Cho_TstScreen, Cho_UseAbs, Damp, DiaMax, DiaMin, DiaMnZ, Did_DecDrv, DSPNm, DSubScr, Frac_ChVBuf, HaltIt, &
                    iABMnZ, iAlQua, iAtomShl, iBas, iBasSh, iChkQ, IFCSew, iiBstR, iiBstRSh, IndRed, IndRSh, InfRed, InfVec, &
                    IntMap, iOff_Col, iOffQ, iQuAB, iRS2F, iScr, iShlSO, iSOShl, iSP2F, lBuf, LuCho, LuMap, LuPri, LuRed, LuRst, &
                    LuScr, LuSel, LuTmp, MaxQual, MaxRed, MaxVec, MinQual, mmBstRT, Mode_Screen, ModRst, Mx2Sh, MxORSh, MxShPr, &
                    N1_Qual, N1_VecRd, N2_Qual, N2_VecRd, N_Subtr, nBas, nBasSh, nBasT, nBstSh, nChkQ, nCol_Chk, nColAB, nDecom, &
                    nDGM_call, nDimRS, nInteg, nnBstR, nnBstRSh, nnBstRT, nnShl, nnShl_Tot, nnZTot, nQual, nSection, nShell, nSym, &
                    nSys_call, NumCho, NumChT, nVecRS1, RstCho, RstDia, Run_Mode, ScDiag, ShA, ShAB, ShB, ShC, ShCD, ShD, Span, &
                    SSTau, SubScrStat, tDecDrv, tDecom, Thr_PreScreen, ThrCom, ThrDef, ThrDiag, ThrNeg, TimSec, tInteg, &
                    Tol_DiaChk, TooNeg, Trace_Idle, WarNeg, XCho_AdrVec, XDamp, XlDiag, XnBas, XnnShl, XnPass, XnShell, XnSym, &
                    XScDiag, XSpan, XThrCom, XThrDiag, XThrNeg, XTooNeg, XWarNeg
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: Lunit
integer(kind=iwp) :: i, j
character(len=*), parameter :: SecNam = 'Cho_Dump'

irc = 0

write(Lunit,*)
write(Lunit,*)
write(Lunit,*) '>>> Output from ',SecNam,':'
write(Lunit,*)
call XFlush(Lunit)

write(Lunit,*) '*** Contents of Choleksy:'
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'iBas : ',(iBas(i),i=1,8)
write(Lunit,*) 'nBas : ',(nBas(i),i=1,8)
write(Lunit,*) 'XnBas: ',(XnBas(i),i=1,8)
write(Lunit,*) 'nBasT: ',nBasT
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'ThrDef        : ',ThrDef
write(Lunit,*) 'ThrCom        : ',ThrCom
write(Lunit,*) 'ThrDiag       : ',ThrDiag
write(Lunit,*) 'Tol_DiaChk    : ',Tol_DiaChk
write(Lunit,*) 'ThrNeg        : ',ThrNeg
write(Lunit,*) 'WarNeg        : ',WarNeg
write(Lunit,*) 'TooNeg        : ',TooNeg
write(Lunit,*) 'nSym          : ',nSym
write(Lunit,*) 'lBuf          : ',lBuf
write(Lunit,*) 'MinQual       : ',MinQual
write(Lunit,*) 'MaxQual       : ',MaxQual
write(Lunit,*) 'IFCSew        : ',IFCSew
write(Lunit,*) 'Mode_Screen   : ',Mode_Screen
write(Lunit,*) 'Cho_DecAlg    : ',Cho_DecAlg
write(Lunit,*) 'Cho_DecAlg_Def: ',Cho_DecAlg_Def
write(Lunit,*) 'iAlQua        : ',iAlQua
write(Lunit,*) 'MxShPr        : ',MxShPr
write(Lunit,*) 'ModRst        : ',ModRst
write(Lunit,*) 'Run_Mode      : ',Run_Mode
write(Lunit,*) 'ScDiag        : ',ScDiag
write(Lunit,*) 'ChkOnly       : ',ChkOnly
write(Lunit,*) 'Cho_IntChk    : ',Cho_IntChk
write(Lunit,*) 'Cho_MinChk    : ',Cho_MinChk
write(Lunit,*) 'Cho_UseAbs    : ',Cho_UseAbs
write(Lunit,*) 'Cho_TrcNeg    : ',Cho_TrcNeg
write(Lunit,*) 'Cho_ReOrd     : ',Cho_ReOrd
write(Lunit,*) 'Cho_DiaChk    : ',Cho_DiaChk
write(Lunit,*) 'Cho_TstScreen : ',Cho_TstScreen
write(Lunit,*) 'Cho_1Center   : ',Cho_1Center
write(Lunit,*) 'Cho_No2Center : ',Cho_No2Center
write(Lunit,*) 'Cho_PreScreen : ',Cho_PreScreen
write(Lunit,*) 'Cho_SimP      : ',Cho_SimP
write(Lunit,*) 'Cho_Fake_Par  : ',Cho_Fake_Par
write(Lunit,*) 'RstDia        : ',RstDia
write(Lunit,*) 'RstCho        : ',RstCho
write(Lunit,*) 'Did_DecDrv    : ',Did_DecDrv
write(Lunit,*) 'HaltIt        : ',HaltIt
write(Lunit,*) 'Trace_Idle    : ',Trace_Idle
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'LuCho: ',(LuCho(i),i=1,8)
write(Lunit,*) 'LuSel: ',(LuSel(i),i=1,8)
write(Lunit,*) 'LuTmp: ',(LuTmp(i),i=1,8)
write(Lunit,*) 'LuPri: ',LuPri
write(Lunit,*) 'LuScr: ',LuScr
write(Lunit,*) 'LuRed: ',LuRed
write(Lunit,*) 'LuRst: ',LuRst
write(Lunit,*) 'LuMap: ',LuMap
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'nShell   : ',nShell
write(Lunit,*) 'nnShl_Tot: ',nnShl_Tot
write(Lunit,*) 'nnShl    : ',nnShl
write(Lunit,*) 'MxORSh   : ',MxORSh
write(Lunit,*) 'Mx2Sh    : ',Mx2Sh
do j=1,3
  write(Lunit,*) 'iiBstR : ',(iiBstR(i,j),i=1,8)
end do
do j=1,3
  write(Lunit,*) 'nnBstR : ',(nnBstR(i,j),i=1,8)
end do
write(Lunit,*) 'nnBstRT: ',(nnBstRT(i),i=1,3)
write(Lunit,*) 'mmBstRT: ',mmBstRT
write(Lunit,*) 'nQual  : ',(nQual(i),i=1,8)
write(Lunit,*) 'iOffQ  : ',(iOffQ(i),i=1,8)
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'DiaMax: ',(DiaMax(i),i=1,8)
write(Lunit,*) 'DiaMin: ',(DiaMin(i),i=1,8)
write(Lunit,*) 'Damp  : ',(Damp(i),i=1,2)
write(Lunit,*) 'Span  : ',Span
write(Lunit,*) 'XlDiag: ',XlDiag
write(Lunit,*) 'DiaMnZ: ',DiaMnZ
write(Lunit,*) 'Thr_PreScreen: ',Thr_PreScreen
write(Lunit,*) 'iABMnZ: ',iABMnZ
write(Lunit,*) 'nnZTot: ',nnZTot
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'NumCho   : ',(NumCho(i),i=1,8)
write(Lunit,*) 'NumChT   : ',NumChT
write(Lunit,*) 'MaxVec   : ',MaxVec
write(Lunit,*) 'MaxRed   : ',MaxRed
write(Lunit,*) 'BlockSize: ',BlockSize
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'ShA     : ',ShA
write(Lunit,*) 'ShB     : ',ShB
write(Lunit,*) 'ShAB    : ',ShAB
write(Lunit,*) 'ShC     : ',ShC
write(Lunit,*) 'ShD     : ',ShD
write(Lunit,*) 'ShCD    : ',ShCD
write(Lunit,*) 'nColAB  : ',nColAB
write(Lunit,*) 'iOff_Col: ',(iOff_Col(i),i=1,8)
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'XThrCom    : ',XThrCom
write(Lunit,*) 'XThrDiag   : ',XThrDiag
write(Lunit,*) 'XDamp      : ',(XDamp(i),i=1,2)
write(Lunit,*) 'XSpan      : ',XSpan
write(Lunit,*) 'XThrNeg    : ',XThrNeg
write(Lunit,*) 'XWarNeg    : ',XWarNeg
write(Lunit,*) 'XTooNeg    : ',XTooNeg
write(Lunit,*) 'XnSym      : ',XnSym
write(Lunit,*) 'XnShell    : ',XnShell
write(Lunit,*) 'XnnShl     : ',XnnShl
write(Lunit,*) 'XnPass     : ',XnPass
write(Lunit,*) 'XCho_AdrVec: ',XCho_AdrVec
write(Lunit,*) 'XScDiag    : ',XScDiag
write(Lunit,*)
call XFlush(Lunit)
do j=1,nChkQ+1
  write(Lunit,*) 'iChkQ   : ',(iChkQ(i,j),i=1,4)
end do
write(Lunit,*) 'nCol_Chk: ',nCol_Chk
do j=1,nSection
  write(Lunit,*) 'TimSec  : ',(TimSec(i,j),i=1,4)
end do
do j=1,nInteg
  write(Lunit,*) 'tInteg  : ',(tInteg(i,j),i=1,2)
end do
do j=1,nDecom
  write(Lunit,*) 'tDecom  : ',(tDecom(i,j),i=1,2)
end do
write(Lunit,*) 'tDecDrv : ',(tDecDrv(i),i=1,2)
write(Lunit,*)
write(Lunit,*) 'nVecRS1: ',(nVecRS1(i),i=1,8)
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'Cho_AdrVec: ',Cho_AdrVec
write(Lunit,*) 'Cho_IOVec : ',Cho_IOVec
write(Lunit,*) 'nSys_Call : ',nSys_Call
write(Lunit,*) 'nDGM_Call : ',nDGM_Call
write(Lunit,*) 'N1_VecRd  : ',N1_VecRd
write(Lunit,*) 'N2_VecRd  : ',N2_VecRd
write(Lunit,*) 'N_Subtr   : ',N_Subtr
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'N1_Qual: ',N1_Qual
write(Lunit,*) 'N2_Qual: ',N2_Qual
write(Lunit,*)
write(Lunit,*) 'Frac_ChVBuf: ',Frac_ChVBuf
call XFlush(Lunit)
write(Lunit,*)
write(Lunit,*) '    (dimension)'
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'InfRed  : ',size(InfRed)
write(Lunit,*) 'InfVec  : ',size(InfVec)
write(Lunit,*) 'IndRed  : ',size(IndRed)
write(Lunit,*) 'IndRSh  : ',size(IndRsh)
write(Lunit,*) 'iScr    : ',size(iScr)
write(Lunit,*) 'iiBstRSh: ',size(iiBstRSh)
write(Lunit,*) 'nnBstRSh: ',size(nnBstRSh)
write(Lunit,*) 'IntMap  : ',size(IntMap)
write(Lunit,*) 'nDimRS  : ',size(nDimRS)
write(Lunit,*) 'iRS2F   : ',size(iRS2F)
write(Lunit,*) 'iSOShl  : ',size(iSOShl)
write(Lunit,*) 'iShlSO  : ',size(iShlSO)
write(Lunit,*) 'iQuab   : ',size(iQuab)
write(Lunit,*) 'iBasSh  : ',size(iBasSh)
write(Lunit,*) 'nBasSh  : ',size(nBasSh)
write(Lunit,*) 'nBstSh  : ',size(nBstSh)
write(Lunit,*) 'iAtomShl: ',size(iAtomShl)
write(Lunit,*) 'iSP2F   : ',size(iSP2F)
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*)
write(Lunit,*) 'Cho_SScreen: ',Cho_SScreen
write(Lunit,*) 'SSTau      : ',SSTau
write(Lunit,*) 'SubScrStat : ',(SubScrStat(i),i=1,2)
write(Lunit,*) 'DSubScr    : ',size(DSubScr)
write(Lunit,*) 'DSPNm      : ',size(DSPNm)

end subroutine Cho_Dump
