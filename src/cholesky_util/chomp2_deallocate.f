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
* Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
*               2010, Jonas Bostrom                                    *
*               2021, Roland Lindh                                     *
************************************************************************
      SubRoutine ChoMP2_deallocate(irc)
      use ChoMP2, only: ChoMP2_allocated
      use ChoMP2, only: iFirst, iFirstS, NumOcc, LnOcc, LnT1am, LiT1am
      use ChoMP2, only: LnMatij, LiMatij, lUnit, NumBatOrb, LnBatOrb
      use ChoMP2, only: LnPQprod, LiPQprod
C
C     Purpose: to deallocate memory of the  Cholesky MP2 program.
C
#include "implicit.fh"
#include "chomp2.fh"
#include "stdalloc.fh"

      irc = 0

      Call ChoMP2g_deallocate(irc)

      If (.NOT.ChoMP2_allocated) Return

      Call mma_deallocate(LiPQprod)
      Call mma_deallocate(LnPQprod)
      Call mma_deallocate(LnBatOrb)
      Call mma_deallocate(NumBatOrb)
      Call mma_deallocate(lUnit)
      Call mma_deallocate(LiMatij)
      Call mma_deallocate(LnMatij)
      Call mma_deallocate(LiT1am)
      Call mma_deallocate(LnT1am)
      Call mma_deallocate(LnOcc)
      Call mma_deallocate(NumOcc)
      Call mma_deallocate(iFirstS)
      Call mma_deallocate(iFirst)
*
      ChoMP2_allocated=.FALSE.

      End

      SubRoutine ChoMP2g_deallocate(irc)
      use ChoMP2, only: ChoMP2g_allocated, EFrozT, EOccuT, EVirtT
      use ChoMP2, only: AdrR1, AdrR2
      use ChoMP2, only: MP2D_full, MP2D
      use ChoMP2, only: MP2W_full, MP2W
      use ChoMP2, only: MP2D_e_full, MP2D_e
      use ChoMP2, only: MP2W_e_full, MP2W_e
*
*     Purpose: Deallocate memory needed for
*              MP2-gradients or properties.
*
#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "stdalloc.fh"

      irc = 0

      If (.NOT.ChoMP2g_allocated) Return

      Call mma_deallocate(MP2D_full)
      Call mma_deallocate(MP2W_full)
      Call mma_deallocate(MP2D_e_full)
      Call mma_deallocate(MP2W_e_full)
      Do iSym = 1, 8
         MP2D(iSym)%A=>Null()
         MP2W(iSym)%A=>Null()
         MP2D_e(iSym)%A=>Null()
         MP2W_e(iSym)%A=>Null()
      End Do
      Call mma_deallocate(AdrR2)
      Call mma_deallocate(AdrR1)
      Call mma_deallocate(EVirtT)
      Call mma_deallocate(EOccuT)
      Call mma_deallocate(EFrozT)

      ChoMP2g_allocated=.FALSE.

      End
