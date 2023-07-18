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
      SubRoutine ChoMP2_Tra(COcc,CVir,Diag,DoDiag)
!
!     Thomas Bondo Pedersen, Dec. 2004.
!
!     Purpose: transform Cholesky vectors to (ai) MO basis.
!
      use stdalloc
      Implicit None
      Real*8  COcc(*), CVir(*), Diag(*)
      Logical DoDiag
#include "cholesky.fh"
#include "chomp2.fh"

      Character(LEN=10), Parameter:: SecNam = 'ChoMP2_Tra'

      Integer kOffD, iSym, lW
      Real*8, Allocatable:: TraMax(:)

!     Allocate remaining memory.
!     --------------------------

      Call mma_maxDBLE(lw)
      Call mma_allocate(TraMax,lw,Label='TraMax')

      kOffD = 1
      Do iSym = 1,nSym

!        Open files for MO vectors.
!        --------------------------

         Call ChoMP2_OpenF(1,1,iSym)

!        Transform vectors.
!        ------------------

         Call ChoMP2_Tra_1(COcc,CVir,Diag(kOffD),DoDiag,TraMax,lW,iSym)
         If (DoDiag) kOffD = kOffD + nT1am(iSym)

!        Close files for MO vectors.
!        ---------------------------

         Call ChoMP2_OpenF(2,1,iSym)

      End Do

!     Free memory.
!     ------------

      Call mma_deallocate(TraMax)

      End
