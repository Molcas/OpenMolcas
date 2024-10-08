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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

subroutine GMFree()

use InfSCF, only: CMO, CMO_ref, Darwin, Dens, EDFT, EOrb, FockAO, FockMO, HDiag, KntE, MssVlc, OccNo, OneHam, OrbType, Ovrlp, TrM, &
                  TwoHam, Vxc
#ifdef _FDE_
use Embedding_Global, only: embInt
#endif
use stdalloc, only: mma_deallocate

implicit none

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Deallocate memory
call mma_deallocate(Darwin)
call mma_deallocate(MssVlc)
call mma_deallocate(KntE)
call mma_deallocate(EDFT)
call mma_deallocate(TwoHam)
call mma_deallocate(Vxc)
call mma_deallocate(Dens)
call mma_deallocate(OrbType)
call mma_deallocate(EOrb)
call mma_deallocate(OccNo)
call mma_deallocate(FockMO)
call mma_deallocate(FockAO)
call mma_deallocate(CMO_ref)
call mma_deallocate(CMO)
call mma_deallocate(TrM)

call mma_deallocate(Ovrlp)
call mma_deallocate(OneHam)
call mma_deallocate(HDiag)
#ifdef _FDE_
call mma_deallocate(embInt,safe='*')
#endif

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine GMFree
