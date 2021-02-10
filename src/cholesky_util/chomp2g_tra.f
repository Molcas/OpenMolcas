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
* Copyright (C) 2004, Jonas Bostrom                                    *
************************************************************************
      SubRoutine ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag,iMoType1,iMoType2)
C
C     Jonas Bostrom, Dec. 2004.
C
C     Purpose: transform Cholesky vectors to (pq) MO basis.
C
#include "implicit.fh"
      Real*8  COrb1(*), COrb2(*), Diag(*)
      Logical DoDiag
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2g.fh"
#include "stdalloc.fh"

      Character(LEN=10), Parameter:: SecNam = 'ChoMP2_Tra'

      Real*8, Allocatable:: TraMax(:)

C     Allocate remaining memory.
C     --------------------------

*     Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
      iVecType = iMoType2 + (iMoType1-1)*nMoType

      Call mma_maxDBLE(lw)
      Call mma_allocate(TraMax,lW,Label='TraMax')

      kOffD = 1
      Do iSym = 1,nSym

C        Open files for MO vectors.
C        --------------------------

         Call ChoMP2_OpenF(1,1,iSym)

C        Transform vectors.
C        ------------------

         Call ChoMP2g_Tra_1(COrb1,COrb2,Diag(kOffD),DoDiag,TraMax,lW,
     &                     iSym,iMoType1,iMoType2)
         kOffD = kOffD + nMoMo(iSym,iVecType)

C        Close files for MO vectors.
C        ---------------------------

         Call ChoMP2_OpenF(2,1,iSym)

      End Do

C     Free memory.
C     ------------

      Call mma_deallocate(TraMax)

      End
