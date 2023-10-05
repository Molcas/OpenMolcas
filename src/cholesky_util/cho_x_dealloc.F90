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

use Cholesky, only: iAtomShl, iBasSh, iiBstRSh, iiBstRSh_G, iiBstRSh_Hidden, iiBstRSh_L_Hidden, iL2G, IndRed, IndRed_G, &
                    IndRed_G_Hidden, IndRed_Hidden, IndRSh, IndRSh_G, IndRSh_G_Hidden, IndRSh_Hidden, InfRed, InfRed_G, &
                    InfRed_G_Hidden, InfRed_Hidden, InfVec, InfVec_Bak, InfVec_G, InfVec_G_Hidden, InfVec_Hidden, IntMap, iQL2G, &
                    iQuAB, iQuAB_Hidden, iQuAB_L, iQuAB_L_Hidden, iRS2F, iScr, iShlSO, iShP2Q, iShP2RS, iSimRI, iSOShl, iSP2F, &
                    LQ_Tot, nBasSh, nBstSh, nDimRS, nnBstRSh, nnBstRSh_G, nnBstRSh_Hidden, nnBstRSh_L_Hidden
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc

irc = 0

! Deallocate.
! -----------

if (allocated(InfRed_Hidden)) call mma_deallocate(InfRed_Hidden)
nullify(InfRed)

if (allocated(InfVec_Hidden)) call mma_deallocate(InfVec_Hidden)
nullify(InfVec)

if (allocated(IndRed_Hidden)) call mma_deallocate(IndRed_Hidden)
nullify(IndRed)

if (allocated(IndRSh_Hidden)) call mma_deallocate(IndRSh_Hidden)
nullify(IndRSh)

if (allocated(iScr)) call mma_deallocate(iScr)

if (allocated(iiBstRSh_Hidden)) call mma_deallocate(iiBstRSh_Hidden)
nullify(iiBstRSh)

if (allocated(nnBstRSh_Hidden)) call mma_deallocate(nnBstRSh_Hidden)
nullify(nnBstRSh)

if (allocated(IntMap)) call mma_deallocate(IntMap)

if (allocated(nDimRS)) call mma_deallocate(nDimRS)

if (allocated(iRS2F)) call mma_deallocate(iRS2F)

if (allocated(iSOShl)) call mma_deallocate(iSOShl)

if (allocated(iShlSO)) call mma_deallocate(iShlSO)

if (allocated(iQuAB_Hidden)) call mma_deallocate(iQuAB_Hidden)
nullify(iQuAB)

if (allocated(iBasSh)) call mma_deallocate(iBasSh)

if (allocated(nBasSh)) call mma_deallocate(nBasSh)

if (allocated(nBstSh)) call mma_deallocate(nBstSh)

if (allocated(iAtomShl)) call mma_deallocate(iAtomShl)

if (allocated(iSP2F)) call mma_deallocate(iSP2F)

if (allocated(iShP2RS)) call mma_deallocate(iShP2RS)

if (allocated(iShP2Q)) call mma_deallocate(iShP2Q)

if (allocated(iQuAB_L_Hidden)) call mma_deallocate(iQuAB_L_Hidden)
nullify(iQuAB_L)

if (allocated(iQL2G)) call mma_deallocate(iQL2G)

if (allocated(LQ_Tot)) call mma_deallocate(LQ_Tot)

if (allocated(InfVec_Bak)) call mma_deallocate(InfVec_Bak)

if (allocated(iSimRI)) call mma_deallocate(iSimRI)

if (allocated(InfVec_G_Hidden)) call mma_deallocate(InfVec_G_Hidden)
nullify(InfVec_G)

if (allocated(IndRed_G_Hidden)) call mma_deallocate(IndRed_G_Hidden)
nullify(IndRed_G)

if (allocated(InfRed_G_Hidden)) call mma_deallocate(InfRed_G_Hidden)
nullify(InfRed_G)

if (allocated(IndRSh_G_Hidden)) call mma_deallocate(IndRSh_G_Hidden)
nullify(IndRSh_G)

if (allocated(iiBstRSh_L_Hidden)) call mma_deallocate(iiBstRSh_L_Hidden)
nullify(iiBstRSh_G)

if (allocated(nnBstRSh_L_Hidden)) call mma_deallocate(nnBstRSh_L_Hidden)
nullify(nnBstRSh_G)

if (allocated(iL2G)) call mma_deallocate(iL2G)

return

end subroutine Cho_X_Dealloc
