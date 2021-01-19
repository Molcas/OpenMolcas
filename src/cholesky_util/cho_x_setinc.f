************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Cho_X_SetInc(irc)
C
C     T.B. Pedersen, July 2004.
C
C     Purpose: define all entries in include files
C              choprint.fh
C              choorb.fh
C              cholesky.fh
C              choptr2.fh
C              chosew.fh
C              cholq.fh
C              chovecbuf.fh
C              chosubscr.fh
C              chosimri.fh
C              chopar.fh
C              cho_para_info.fh
C              chobkm.fh
C
      Implicit None
      Integer irc
#include "choorb.fh"
#include "choprint.fh"
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "choptr2.fh"
#include "chosew.fh"
#include "cholq.fh"
#include "chosubscr.fh"
#include "chosimri.fh"
#include "chopar.fh"
#include "cho_para_info.fh"
#include "chobkm.fh"

      Integer iLarge
      Parameter (iLarge = 99999999)

      Real*8 Large, Small
      Parameter (Large = 1.0D15, Small = 1.0D-15)

C     Set return code.
C     ----------------

      irc = 0

C     choprint.fh.
C     -------------

      iPrint = -iLarge

C     choorb.fh.
C     -----------

      Call Cho_iZero(iBas,8)
      Call Cho_iZero(nBas,8)
      Call Cho_iZero(XnBas,8)
      nBasT = 0

C     cholesky.fh.
C     -------------

      ThrCom  = Large
      ThrDiag = Large
      Tol_DiaChk = -Large
      ThrNeg  = Large
      WarNeg  = Large
      TooNeg  = Large
      nSym    = -iLarge
      lBuf    = -iLarge
      MinQual = iLarge
      MaxQUal = iLarge
      IFCSew  = -iLarge
      Mode_Screen = -iLarge
      Cho_DecAlg  = -iLarge
      Cho_DecAlg_Def = -iLarge
      iAlQua  = iLarge
      MxShPr  = iLarge
      ModRst  = iLarge
      Run_Mode = iLarge
      ScDiag  = .false.
      ChkOnly = .false.
      Cho_IntChk = .false.
      Cho_MinChk = .false.
      Cho_UseAbs = .false.
      Cho_TrcNeg = .false.
      Cho_ReOrd  = .false.
      Cho_DiaChk = .false.
      Cho_TstScreen = .false.
      Cho_1Center = .false.
      Cho_No2Center = .false.
      Cho_PreScreen = .false.
      Cho_SimP = .false.
      Cho_Fake_Par = .false.
      RstDia  = .false.
      RstCho  = .false.
      Did_DecDrv = .false.
      HaltIt  = .false.
      Trace_Idle = .false.

      Call Cho_iZero(LuCho,8)
      Call Cho_iZero(LuSel,8)
      Call Cho_iZero(LuTmp,8)
      LuPri = 0
      LuScr = 0
      LuRed = 0
      LuRst = 0
      LuMap = 0

      nShell = 0
      nnShl_Tot  = 0
      nnShl = 0
      MxORSh = 0
      Mx2Sh  = 0
      Call Cho_iZero(iiBstR,8*3)
      Call Cho_iZero(nnBstR,8*3)
      Call Cho_iZero(nnBstRT,3)
      mmBstRT = 0
      Call Cho_iZero(nQual,8)
      Call Cho_iZero(iOffQ,8)

      Call Cho_dZero(DiaMax,8)
      Call Cho_dZero(DiaMaxT,8)
      Call Cho_dZero(DiaMin,8)
      Call Cho_dZero(Damp,2)
      Span   = Large
      XlDiag = Large
      DiaMnZ = Large
      Thr_PreScreen = -Large
      iABMnZ = -iLarge
      nnZTot = 0

      Call Cho_iZero(NumCho,8)
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
      Call Cho_iZero(iOff_Col,8)

      XThrCom  = Large
      XThrDiag = Large
      Call Cho_dZero(XDamp,2)
      XSpan    = Large
      XThrNeg  = Large
      XWarNeg  = Large
      XTooNeg  = Large
      XnSym    = 0
      XnShell  = 0
      XnnShl   = 0
      XnPass   = 0
      XScDiag  = .false.
      XCho_AdrVec = -iLarge

      Call Cho_iZero(iChkQ,4*(nChkQ+1))
      nCol_Chk = -iLarge
      Call Cho_dZero(TimSec,4*nSection)
      Call Cho_dZero(tInteg,2*nInteg)
      Call Cho_dZero(tDecom,2*nDecom)
      Call Cho_dZero(tMisc,2*nMisc)
      Call Cho_dZero(tDecDrv,2)

      Call Cho_iZero(nVecRS1,8)

      Cho_AdrVec= -iLarge
      Cho_IOVec = -iLarge
      nSys_Call = 0
      nDGM_Call = 0
      N1_VecRd  = 0
      N2_VecRd  = 0
      N_Subtr   = 0

      N1_Qual = -iLarge
      N2_Qual = iLarge

      Frac_ChVBuf = 0.0d0

C     Zero memory in pointers in chosew.fh.
C     --------------------------------------

      ip_iShP2RS = 0
      l_iShP2RS  = 0
      ip_iShP2Q = 0
      l_iShP2Q  = 0
      ip_iOff_Batch = 0
      l_iOff_Batch  = 0
      Call Cho_iZero(nDim_Batch,8)

C     cholq.fh.
C     ----------

      Call Cho_iZero(nQual_L,8)
      Call Cho_iZero(ip_LQ_Sym,8)
      Call Cho_iZero(l_LQ_Sym,8)
      Call Cho_iZero(ldLQ,8)
      ip_iQL2G = 0
      l_iQL2G  = 0
      ip_LQ = 0
      l_LQ  = 0

C     Zero memory pointers in choptr2.fh.
C     ------------------------------------

      ip_mySP=0
      l_mySP=0
      n_mySP=0
      ip_Idle=0
      l_Idle=0

C     chovecbuf.fh.
C     --------------

      ip_ChVBuf = 0
      l_ChvBuf  = 0
      ip_ChVBfI = 0
      l_ChvBfI  = 0
      Call Cho_iZero(ip_ChVBuf_Sym,8)
      Call Cho_iZero(l_ChVBuf_Sym,8)
      Call Cho_iZero(ip_ChVBfI_Sym,8)
      Call Cho_iZero(l_ChVBfI_Sym,8)
      Call Cho_iZero(nVec_in_Buf,8)

C     chosubscr.fh.
C     --------------

      Cho_SScreen = .false.
      SSTau       = 0.0d0
      SubScrStat(1) = 0.0d0
      SubScrStat(2) = 0.0d0
      ip_DSubScr  = 0
      l_DSubScr   = 0
      ip_DSPNm    = 0
      l_DSPNm     = 0
      SSNorm      = 'tbp'

C     chosimri.fh.
C     -------------

      Cho_SimRI = .false.
      ip_iSimRI = 0
      l_iSimRI  = 0
      Thr_SimRI = -Large

C     chopar.fh.
C     -----------

      Call Cho_iZero(NumCho_Bak,8)

C     cho_para_info.fh.
C     ------------------

      Cho_Real_Par = .false.

C     chobkm.fh.
C     -----------

      ip_BkmVec=0
      l_BkmVec=0
      nRow_BkmVec=0
      nCol_BkmVec=0
      ip_BkmThr=0
      l_BkmThr=0
      nRow_BkmThr=0
      nCol_BkmThr=0

      End
