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

! saved info for parallel runs (chpari):
!
! NumCho_Bak

use Data_Structures, only: V2
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ip_CHVBFI_SYM(8), ip_CHVBUF_SYM(8), l_CHVBFI_SYM(8), l_CHVBUF_SYM(8), LUFV(8,8), n_MySP, NABPK(8,8), &
                     nCol_BkmThr = 0, nCol_BkmVec = 0, nDim_Batch(8), NNBST(8), nnShl_SP, nQual_L(8), nRow_BkmThr = 0, &
                     nRow_BkmVec = 0, NumCho_Bak(8), nVec_in_Buf(8)
real(kind=wp) :: SSTau, SubScrStat(2)
logical(kind=iwp) :: Cho_SScreen
character(len=3) :: SSNorm
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
real(kind=wp), allocatable, target :: BkmThr(:,:), CHVBFI(:), CHVBUF(:), Diag_G_Hidden(:), Diag_Hidden(:), LQ_Tot(:)
real(kind=wp), pointer :: Diag(:) => null(), Diag_G(:) => null()
integer(kind=iwp), parameter :: ChoIniCheck = -6543210
character(len=*), parameter :: REONAM = 'CHFV'

public :: BkmThr, BkmVec, Cho_SScreen, ChoIniCheck, CHVBFI, CHVBUF, Diag, Diag_G, Diag_G_Hidden, Diag_Hidden, DSPNm, DSubScr, &
          iAtomShl, iBasSh, Idle, iiBstRSh, iiBstRSh_G, iiBstRSh_Hidden, iiBstRSh_L_Hidden, iL2G, IndRed, IndRed_G, &
          IndRed_G_Hidden, IndRed_Hidden, IndRSh, IndRSh_G, IndRSh_G_Hidden, IndRSh_Hidden, InfRed, InfRed_G, InfRed_G_Hidden, &
          InfRed_Hidden, InfVec, InfVec_Bak, InfVec_G, InfVec_G_Hidden, InfVec_Hidden, IntMap, iOff_Batch, ip_CHVBFI_SYM, &
          ip_CHVBUF_SYM, iQL2G, iQuAB, iQuAB_here, iQuAB_Hidden, iQuAB_L, iQuAB_L_Hidden, iRS2F, iScr, iShlSO, iShP2Q, iShP2RS, &
          iSimRI, iSOShl, iSP2F, l_CHVBFI_SYM, l_CHVBUF_SYM, LQ, LQ_Tot, LUFV, MySP, n_MySP, NABPK, nBasSh, nBstSh, nCol_BkmThr, &
          nCol_BkmVec, nDim_Batch, nDimRS, NNBST, nnBstRSh, nnBstRSh_G, nnBstRSh_Hidden, nnBstRSh_L_Hidden, nnShl_SP, nQual_L, &
          nRow_BkmThr, nRow_BkmVec, NumCho_Bak, nVec_in_Buf, pTemp, pTemp1, pTemp3, REONAM, SSNorm, SSTau, SubScrStat

end module Cholesky
