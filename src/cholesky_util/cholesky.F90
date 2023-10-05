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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module Cholesky

! For Cholesky decomposition program.
!
! BLOCKSIZE, CHKONLY, CHO_1CENTER, CHO_ADRVEC, CHO_DECALG, CHO_DECALG_DEF, CHO_DIACHK, CHO_FAKE_PAR, CHO_INTCHK, CHO_IOVEC,
! CHO_MINCHK, CHO_NDECALG, CHO_NO2CENTER, CHO_PRESCREEN, CHO_REORD, CHO_SIMP, CHO_TRCNEG, CHO_TSTSCREEN, CHO_USEABS, DAMP, DIAMAX,
! DIAMAXT DIAMIN, DIAMNZ, DID_DECDRV, FRAC_CHVBUF, HALTIT, IABMNZ, IALQUA, ICHKQ, IFCSEW, IIBSTR, INFVEC_N2, IOFF_COL, IOFFQ, LBUF,
! LUCHO, LUMAP LUPRI, LURED, LURST, LUSCR, LUSEL, LUTMP, MAXQUAL, MAXRED, MAXVEC, MINQUAL, MMBSTRT, MODE_SCREEN, MODRST, MX2SH,
! MXORSH, MXSHPR N1_QUAL, N1_VECRD, N2_QUAL, N2_VECRD, N_SUBTR, NALQUA, NCHKQ, NCOL_CHK, NCOLAB, NDECOM, NDGM_CALL, NINTEG, NMISC,
! NNBSTR, NNBSTRT, NNSHL, NNSHL_TOT, NNZTOT, NQUAL NSECTION, NSHELL, NSYM, NSYS_CALL, NUMCHO, NUMCHT, NVECRS1, RSTCHO, RSTDIA,
! RUN_EXTERNAL, RUN_INTERNAL, RUN_MODE, SCDIAG, SHA, SHAB, SHB SHC, SHCD, SHD, SPAN, TDECDRV, TDECOM, THR_PRESCREEN, THRCOM,
! THRDEF, THRDIAG, THRNEG, TIMSEC, TINTEG, TMISC, TOL_DIACHK, TOONEG TRACE_IDLE, WARNEG, XCHO_ADRVEC, XDAMP, XLDIAG, XNNSHL,
! XNPASS, XNSHELL, XNSYM, XSCDIAG, XSPAN, XTHRCOM, XTHRDIAG, XTHRNEG XTOONEG, XWARNEG

! Additional stuff needed for simulation of RI.
!
! CHO_SIMRI, THR_SIMRI

! ChoIniCheck is ONLY for internal use by the Cholesky (external)
! initialization and finalization routines (choini)

! Stuff for Cholesky vector subtraction screening (chosubscr):
!
! Cho_SScreen, DSPNm, DSubScr, SSNorm, SSTau, SubScrStat

! info for Cholesky bookmarks (chobkm):
!
! BkmThr, BkmVec, nCol_BkmThr, nCol_BkmVec, nRow_BkmThr, nRow_BkmVec

! Data for Cholesky vector reordering (choreo):
!
! LUFV, NABPK, NNBST, REONAM

! Information about the Cholesky vector buffer (chovecbuf):
!
! CHVBFI, CHVBUF, ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, l_CHVBUF_SYM, nVec_in_Buf

! Cholesky orbital information
!
! IBAS, NBAS, NBAST, XNBAS

! ChFracMem: fraction of memory used as a buffer for keeping Cholesky vectors previously read from disk (chopar)

! Information for printing during decomposition (choprint):
!
! IPRINT       : print level (from input)
! INF_TIMING   : timing info (cho_prtim)
! INF_INIT     : initialization info
! INF_DIAG     : initial diagonal screening info
! INF_VECBUF   : vector buffer info
! INF_PASS     : pass info
! INF_IN2      : col. info (summary)
! INF_INT      : shell quadr. info
! INF_PROGRESS : decom. tables for each pass
! INF_SUBTR1   : timing info in cho_subtr1 (debug)

! Global stuff for parallel Cholesky
!
! iiBstR_G, LuCho_G, LuRed_G, LuRst_G, mmBstRT_G, myNumCho, nLoc_G, nnBstR_G, nnBstRT_G, nnShl_G, NumCho_G, NumChT_G

! saved info for parallel runs (chpari):
!
! NumCho_Bak

! Cholesky parallel info: additional information if running FAKE parallel,
! which seems to be a special Cholesky parallel mode for e.g. debugging:
!
! Cho_Real_Par

use Data_Structures, only: V2
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: NCHKQ = 12, NDECOM = 4, NINTEG = 2, nLoc_G = 3, NMISC = 5, NSECTION = 9
integer(kind=iwp) :: BLOCKSIZE, CHO_ADRVEC, CHO_DECALG, CHO_DECALG_DEF, CHO_IOVEC, IABMNZ, IALQUA, iBas(8), ICHKQ(4,NCHKQ+1), &
                     IFCSEW, IIBSTR(8,3), iiBstR_G(8,nLoc_G), IOFF_COL(8), IOFFQ(8), ip_CHVBFI_SYM(8), ip_CHVBUF_SYM(8), IPRINT, &
                     l_CHVBFI_SYM(8), l_CHVBUF_SYM(8), LBUF, LUCHO(8), LuCho_G(8), LUFV(8,8), LUMAP, LUPRI, LURED, LuRed_G, LURST, &
                     LuRst_G, LUSCR, LUSEL(8), LUTMP(8), MAXQUAL, MAXRED, MAXVEC, MINQUAL, MMBSTRT, mmBstRT_G, MODE_SCREEN, &
                     MODRST, MX2SH, MXORSH, MXSHPR, myNumCho(8), N1_QUAL, N1_VECRD, N2_QUAL, N2_VECRD, n_MySP, N_SUBTR, &
                     NABPK(8,8), nBas(8), nBasT, nCol_BkmThr = 0, nCol_BkmVec = 0, NCOL_CHK, NCOLAB, NDGM_CALL, nDim_Batch(8), &
                     NNBST(8), NNBSTR(8,3), nnBstR_G(8,nLoc_G), NNBSTRT(3), nnBstRT_G(nLoc_G), NNSHL, nnShl_G, nnShl_SP, &
                     NNSHL_TOT, NNZTOT, NQUAL(8), nQual_L(8), nRow_BkmThr = 0, nRow_BkmVec = 0, NSHELL, NSYM, NSYS_CALL, &
                     NUMCHO(8), NumCho_Bak(8), NumCho_G(8), NUMCHT, NumChT_G, nVec_in_Buf(8), NVECRS1(8), RUN_MODE, SHA, SHAB, &
                     SHB, SHC, SHCD, SHD, XCHO_ADRVEC, XnBas(8), XNNSHL, XNPASS, XNSHELL, XNSYM
real(kind=wp) :: ChFracMem, DAMP(2), DIAMAX(8), DIAMAXT(8), DIAMIN(8), DIAMNZ, FRAC_CHVBUF, SPAN, SSTau, SubScrStat(2), &
                 TDECDRV(2), TDECOM(2,NDECOM), THR_PRESCREEN, THR_SIMRI, THRCOM, THRDIAG, THRNEG, TIMSEC(4,NSECTION), &
                 TINTEG(2,NINTEG), TMISC(2,NMISC), TOL_DIACHK, TOONEG, WARNEG, XDAMP(2), XLDIAG, XSPAN, XTHRCOM, XTHRDIAG, &
                 XTHRNEG, XTOONEG, XWARNEG
logical(kind=iwp) :: CHKONLY, CHO_1CENTER, CHO_DIACHK, CHO_FAKE_PAR, CHO_INTCHK, CHO_MINCHK, CHO_NO2CENTER, CHO_PRESCREEN, &
                     Cho_Real_Par, CHO_REORD, CHO_SIMP, CHO_SIMRI, Cho_SScreen, CHO_TRCNEG, CHO_TSTSCREEN, CHO_USEABS, DID_DECDRV, &
                     HALTIT, RSTCHO, RSTDIA, SCDIAG, timings, TRACE_IDLE, XSCDIAG
character(len=3) :: SSNorm, tv2disk
type(V2) :: LQ(8)
integer(kind=iwp), allocatable :: BkmVec(:,:), iAtomShl(:), iBasSh(:,:), Idle(:), iL2G(:), IntMap(:), iOff_Batch(:,:), iQL2G(:,:), &
                                  iRS2F(:,:), iScr(:), iShlSO(:), iShP2Q(:,:), iShP2RS(:,:), iSimRI(:), iSOShl(:), iSP2F(:), &
                                  MySP(:), nBasSh(:,:), nBstSh(:), nDimRS(:,:)
integer(kind=iwp), allocatable, target :: iiBstRSh_Hidden(:,:,:), iiBstRSh_L_Hidden(:,:,:), IndRed_G_Hidden(:,:), &
                                          IndRed_Hidden(:,:), IndRSh_G_Hidden(:), IndRSh_Hidden(:), InfRed_G_Hidden(:), &
                                          InfRed_Hidden(:), InfVec_Bak(:,:,:), InfVec_G_Hidden(:,:,:), InfVec_Hidden(:,:,:), &
                                          iQuAB_here(:,:), iQuAB_Hidden(:,:), iQuAB_L_Hidden(:,:), nnBstRSh_Hidden(:,:,:), &
                                          nnBstRSh_L_Hidden(:,:,:)
integer(kind=iwp), pointer :: iiBstRSh(:,:,:) => null(), iiBstRSh_G(:,:,:) => null(), IndRed(:,:) => null(), &
                              IndRed_G(:,:) => null(), IndRSh(:) => null(), IndRSh_G(:) => null(), InfRed(:) => null(), &
                              InfRed_G(:) => null(), InfVec(:,:,:) => null(), InfVec_G(:,:,:) => null(), iQuAB(:,:) => null(), &
                              iQuAB_L(:,:) => null(), nnBstRSh(:,:,:) => null(), nnBstRSh_G(:,:,:) => null(), &
                              pTemp(:,:) => null(), pTemp1(:) => null(), pTemp3(:,:,:) => null()
real(kind=wp), allocatable :: DSPNm(:), DSubScr(:)
real(kind=wp), allocatable, target :: BkmThr(:,:), CHVBFI(:,:), CHVBUF(:), Diag_G_Hidden(:), Diag_Hidden(:), LQ_Tot(:)
real(kind=wp), pointer :: Diag(:) => null(), Diag_G(:) => null()
integer(kind=iwp), parameter :: CHO_NDECALG = 6, ChoIniCheck = -6543210, INF_DIAG = 6, INF_IN2 = 5, INF_INIT = 3, INF_INT = 6, &
                                INF_PASS = 3, INF_PROGRESS = 4, INF_SUBTR1 = 6, INF_TIMING = 2, INF_VECBUF = 3, INFVEC_N2 = 5, &
                                NALQUA = 2, RUN_EXTERNAL = 2, RUN_INTERNAL = 1
real(kind=wp), parameter :: THRDEF = 1.0D-4
character(len=*), parameter :: REONAM = 'CHFV'

public :: BkmThr, BkmVec, BLOCKSIZE, ChFracMem, CHKONLY, CHO_1CENTER, CHO_ADRVEC, CHO_DECALG, CHO_DECALG_DEF, CHO_DIACHK, &
          CHO_FAKE_PAR, CHO_INTCHK, CHO_IOVEC, CHO_MINCHK, CHO_NDECALG, CHO_NO2CENTER, CHO_PRESCREEN, Cho_Real_Par, CHO_REORD, &
          CHO_SIMP, CHO_SIMRI, Cho_SScreen, CHO_TRCNEG, CHO_TSTSCREEN, CHO_USEABS, ChoIniCheck, CHVBFI, CHVBUF, DAMP, Diag, &
          Diag_G, Diag_G_Hidden, Diag_Hidden, DIAMAX, DIAMAXT, DIAMIN, DIAMNZ, DID_DECDRV, DSPNm, DSubScr, FRAC_CHVBUF, HALTIT, &
          IABMNZ, IALQUA, iAtomShl, iBas, iBasSh, ICHKQ, Idle, IFCSEW, IIBSTR, iiBstR_G, iiBstRSh, iiBstRSh_G, iiBstRSh_Hidden, &
          iiBstRSh_L_Hidden, iL2G, IndRed, IndRed_G, IndRed_G_Hidden, IndRed_Hidden, IndRSh, IndRSh_G, IndRSh_G_Hidden, &
          IndRSh_Hidden, INF_DIAG, INF_IN2, INF_INIT, INF_INT, INF_PASS, INF_PROGRESS, INF_SUBTR1, INF_TIMING, INF_VECBUF, InfRed, &
          InfRed_G, InfRed_G_Hidden, InfRed_Hidden, InfVec, InfVec_Bak, InfVec_G, InfVec_G_Hidden, InfVec_Hidden, INFVEC_N2, &
          IntMap, iOff_Batch, IOFF_COL, IOFFQ, ip_CHVBFI_SYM, ip_CHVBUF_SYM, IPRINT, iQL2G, iQuAB, iQuAB_here, iQuAB_Hidden, &
          iQuAB_L, iQuAB_L_Hidden, iRS2F, iScr, iShlSO, iShP2Q, iShP2RS, iSimRI, iSOShl, iSP2F, l_CHVBFI_SYM, l_CHVBUF_SYM, LBUF, &
          LQ, LQ_Tot, LUCHO, LuCho_G, LUFV, LUMAP, LUPRI, LURED, LuRed_G, LURST, LuRst_G, LUSCR, LUSEL, LUTMP, MAXQUAL, MAXRED, &
          MAXVEC, MINQUAL, MMBSTRT, mmBstRT_G, MODE_SCREEN, MODRST, MX2SH, MXORSH, MXSHPR, myNumCho, MySP, N1_QUAL, N1_VECRD, &
          N2_QUAL, N2_VECRD, n_MySP, N_SUBTR, NABPK, NALQUA, nBas, nBasSh, nBasT, nBstSh, NCHKQ, nCol_BkmThr, nCol_BkmVec, &
          NCOL_CHK, NCOLAB, NDECOM, NDGM_CALL, nDim_Batch, nDimRS, NINTEG, nLoc_G, NMISC, NNBST, NNBSTR, nnBstR_G, nnBstRSh, &
          nnBstRSh_G, nnBstRSh_Hidden, nnBstRSh_L_Hidden, NNBSTRT, nnBstRT_G, NNSHL, nnShl_G, nnShl_SP, NNSHL_TOT, NNZTOT, NQUAL, &
          nQual_L, nRow_BkmThr, nRow_BkmVec, NSECTION, NSHELL, NSYM, NSYS_CALL, NUMCHO, NumCho_Bak, NumCho_G, NUMCHT, NumChT_G, &
          nVec_in_Buf, NVECRS1, pTemp, pTemp1, pTemp3, REONAM, RSTCHO, RSTDIA, RUN_EXTERNAL, RUN_INTERNAL, RUN_MODE, SCDIAG, SHA, &
          SHAB, SHB, SHC, SHCD, SHD, SPAN, SSNorm, SSTau, SubScrStat, TDECDRV, TDECOM, THR_PRESCREEN, THR_SIMRI, THRCOM, THRDEF, &
          THRDIAG, THRNEG, timings, TIMSEC, TINTEG, TMISC, TOL_DIACHK, TOONEG, TRACE_IDLE, tv2disk, WARNEG, XCHO_ADRVEC, XDAMP, &
          XLDIAG, XnBas, XNNSHL, XNPASS, XNSHELL, XNSYM, XSCDIAG, XSPAN, XTHRCOM, XTHRDIAG, XTHRNEG, XTOONEG, XWARNEG

end module Cholesky
