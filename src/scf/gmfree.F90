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
      SubRoutine GMFree()
      use SCF_Arrays, only: Darwin, MssVlc, KntE, EDFT, TwoHam, Vxc, Dens, EOrb, OccNo, FockMO, FockAO,     &
                            CMO_ref, CMO, TrM, Ovrlp, OneHam, HDiag
      use Orb_Type, only: OrbType
#ifdef _FDE_
      use Embedding_Global, only: embInt
#endif
      use stdalloc, Only: mma_deallocate
      Implicit None
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!---- Deallocate memory
      Call mma_deallocate(Darwin)
      Call mma_deallocate(MssVlc)
      Call mma_deallocate(KntE)
      Call mma_deallocate(EDFT)
      Call mma_deallocate(TwoHam)
      Call mma_deallocate(Vxc)
      Call mma_deallocate(Dens)
      Call mma_deallocate(OrbType)
      Call mma_deallocate(EOrb)
      Call mma_deallocate(OccNo)
      Call mma_deallocate(FockMO)
      Call mma_deallocate(FockAO)
      Call mma_deallocate(CMO_ref)
      Call mma_deallocate(CMO)
      Call mma_deallocate(TrM)
!
      Call mma_deallocate(Ovrlp)
      Call mma_deallocate(OneHam)
      Call mma_deallocate(HDiag)
#ifdef _FDE_
      If (Allocated(embInt)) Call mma_deallocate(embInt)
#endif
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine GMFree
