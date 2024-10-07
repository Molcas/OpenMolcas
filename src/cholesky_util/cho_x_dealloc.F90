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

call mma_deallocate(InfRed_Hidden,safe='*')

call mma_deallocate(InfVec_Hidden,safe='*')

call mma_deallocate(IndRed_Hidden,safe='*')

call mma_deallocate(IndRSh_Hidden,safe='*')

call mma_deallocate(iScr,safe='*')

call mma_deallocate(iiBstRSh_Hidden,safe='*')

call mma_deallocate(nnBstRSh_Hidden,safe='*')

call mma_deallocate(IntMap,safe='*')

call mma_deallocate(nDimRS,safe='*')

call mma_deallocate(iRS2F,safe='*')

call mma_deallocate(iSOShl,safe='*')

call mma_deallocate(iShlSO,safe='*')

call mma_deallocate(iQuAB_Hidden,safe='*')

call mma_deallocate(iBasSh,safe='*')

call mma_deallocate(nBasSh,safe='*')

call mma_deallocate(nBstSh,safe='*')

call mma_deallocate(iAtomShl,safe='*')

call mma_deallocate(iSP2F,safe='*')

call mma_deallocate(iShP2RS,safe='*')

call mma_deallocate(iShP2Q,safe='*')

call mma_deallocate(iQuAB_L_Hidden,safe='*')

call mma_deallocate(iQL2G,safe='*')

call mma_deallocate(LQ_Tot,safe='*')

call mma_deallocate(InfVec_Bak,safe='*')

call mma_deallocate(iSimRI,safe='*')

call mma_deallocate(InfVec_G_Hidden,safe='*')

call mma_deallocate(IndRed_G_Hidden,safe='*')

call mma_deallocate(InfRed_G_Hidden,safe='*')

call mma_deallocate(IndRSh_G_Hidden,safe='*')

call mma_deallocate(iiBstRSh_L_Hidden,safe='*')

call mma_deallocate(nnBstRSh_L_Hidden,safe='*')

call mma_deallocate(iL2G,safe='*')

nullify(InfRed,InfVec,IndRed,IndRSh,iiBstRSh,nnBstRSh,iQuAB,iQuAB_L,InfVec_G,IndRed_G,InfRed_G,IndRSh_G,iiBstRSh_G,nnBstRSh_G)

return

end subroutine Cho_X_Dealloc
