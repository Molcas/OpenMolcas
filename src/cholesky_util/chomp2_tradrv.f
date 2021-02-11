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
      SubRoutine ChoMP2_TraDrv(irc,CMO,Diag,DoDiag)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: AO-to-MO (ai) transformation of Cholesky vectors
C              performed directly in reduced sets. This assumes
C              that the MP2 program has been appropriately initialized.
C
      Implicit None
      Integer irc
      Real*8  CMO(*), Diag(*)
      Logical DoDiag
#include "cholesky.fh"
#include "chomp2.fh"
#include "stdalloc.fh"

      Character(LEN=6), Parameter:: ThisNm = 'TraDrv'
      Character(LEN=13), Parameter:: SecNam = 'ChoMP2_TraDrv'

      Real*8, Allocatable:: COcc(:), CVir(:)

      irc = 0

C     Reorder MO coefficients.
C     ------------------------

      Call mma_allocate(COcc,nT1AOT(1),Label='COcc')
      Call mma_allocate(CVir,nAOVir(1),Label='CVir')
      Call ChoMP2_MOReOrd(CMO,COcc,CVir)

C     Transform vectors.
C     ------------------

      Call ChoMP2_Tra(COcc,CVir,Diag,DoDiag)

C     Deallocate reordered MO coefficients.
C     -------------------------------------

      Call mma_deallocate(CVir)
      Call mma_deallocate(COcc)

      End
