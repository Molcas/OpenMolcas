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

subroutine Cho_X_Dealloc(irc)
!
! T.B. Pedersen, July 2004.
!
! Purpose: deallocate ALL index arrays for the Cholesky utility.
!          On exit, irc=0 signals sucessful completion.

use Definitions, only: iwp
use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iAtomShl, iShlSO, iRS2F, IntMap, iScr, nDimRS, iL2G, iShP2RS, iShP2Q, &
                  iQL2G, LQ_Tot, iSimRI
use ChoSwp, only: iQuAB, iQuAB_L, iQuAB_Hidden, iQuAB_L_Hidden, nnBstRSh_Hidden, nnBstRSh, nnBstRSh_L_Hidden, nnBstRSh_G, &
                  iiBstRSh_Hidden, iiBstRSh, iiBstRSh_L_Hidden, iiBstRSh_G, IndRSh_Hidden, IndRSh, IndRSh_G_Hidden, IndRSh_G, &
                  InfRed_Hidden, InfRed, InfRed_G_Hidden, InfRed_G, InfVec_Hidden, InfVec, InfVec_G_Hidden, InfVec_G, &
                  IndRed_Hidden, IndRed, IndRed_G_Hidden, IndRed_G, InfVec_Bak
use stdalloc, only: mma_deallocate
use ChPari

implicit none
integer(kind=iwp) :: irc
character(len=13), parameter :: SecNam = 'Cho_X_Dealloc'

irc = 0

! Deallocate.
! -----------

if (allocated(InfRed_Hidden)) call mma_deallocate(InfRed_Hidden)
if (associated(InfRed)) InfRed => null()

if (allocated(InfVec_Hidden)) call mma_deallocate(InfVec_Hidden)
if (associated(InfVec)) InfVec => null()

if (allocated(IndRed_Hidden)) call mma_deallocate(IndRed_Hidden)
if (associated(IndRed)) IndRed => null()

if (allocated(IndRSh_Hidden)) call mma_deallocate(IndRSh_Hidden)
if (associated(IndRSh)) IndRSh => null()

if (allocated(iScr)) call mma_deallocate(iScr)

if (allocated(iiBstRSh_Hidden)) call mma_deallocate(iiBstRSh_Hidden)
if (associated(iiBstRSh)) iiBstRSh => null()

if (allocated(nnBstRSh_Hidden)) call mma_deallocate(nnBstRSh_Hidden)
if (associated(nnBstRSh)) nnBstRSh => null()

if (allocated(IntMap)) call mma_deallocate(IntMap)

if (allocated(nDimRS)) call mma_deallocate(nDimRS)

if (allocated(iRS2F)) call mma_deallocate(iRS2F)

if (allocated(iSOShl)) call mma_deallocate(iSOShl)

if (allocated(iShlSO)) call mma_deallocate(iShlSO)

if (allocated(iQuAB_Hidden)) call mma_deallocate(iQuAB_Hidden)
if (associated(iQuAB)) iQuAB => null()

if (allocated(iBasSh)) call mma_deallocate(iBasSh)

if (allocated(nBasSh)) call mma_deallocate(nBasSh)

if (allocated(nBstSh)) call mma_deallocate(nBstSh)

if (allocated(iAtomShl)) call mma_deallocate(iAtomShl)

if (allocated(iSP2F)) call mma_deallocate(iSP2F)

if (allocated(iShP2RS)) call mma_deallocate(iShP2RS)

if (allocated(iShP2Q)) call mma_deallocate(iShP2Q)

if (allocated(iQuAB_L_Hidden)) call mma_deallocate(iQuAB_L_Hidden)
if (associated(iQuAB_L)) iQuAB_L => null()

if (allocated(iQL2G)) call mma_deallocate(iQL2G)

if (allocated(LQ_Tot)) call mma_deallocate(LQ_Tot)

if (allocated(InfVec_Bak)) call mma_deallocate(InfVec_Bak)

if (allocated(iSimRI)) call mma_deallocate(iSimRI)

if (allocated(InfVec_G_Hidden)) call mma_deallocate(InfVec_G_Hidden)
if (associated(InfVec_G)) InfVec_G => null()

if (allocated(IndRed_G_Hidden)) call mma_deallocate(IndRed_G_Hidden)
if (associated(IndRed_G)) IndRed_G => null()

if (allocated(InfRed_G_Hidden)) call mma_deallocate(InfRed_G_Hidden)
if (associated(InfRed_G)) InfRed_G => null()

if (allocated(IndRSh_G_Hidden)) call mma_deallocate(IndRSh_G_Hidden)
if (associated(IndRSh_G)) IndRSh_G => null()

if (allocated(iiBstRSh_L_Hidden)) call mma_deallocate(iiBstRSh_L_Hidden)
if (associated(iiBstRSh_G)) iiBstRSh_G => null()

if (allocated(nnBstRSh_L_Hidden)) call mma_deallocate(nnBstRSh_L_Hidden)
if (associated(nnBstRSh_G)) nnBstRSh_G => null()

if (allocated(iL2G)) call mma_deallocate(iL2G)

return

end subroutine Cho_X_Dealloc
