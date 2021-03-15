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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************
      SubRoutine ReIndexFrag(Array, nDens, nDens_Valence, nBas,         &
     &                       nBas_Valence,nIrrep)
!***********************************************************************
!                                                                      *
! Input: Array(nDens) filled up to Array(nDens_Valence)                *
! Output: Reindexed Array so fragment densities can be put at the      *
!         proper places.                                               *
! Only needed in case of symmetry.                                     *
!                                                                      *
!***********************************************************************
      Implicit None
#include "real.fh"
      Integer nDens, nDens_Valence,nIrrep
      Real*8  Array(nDens)
      Integer nBas(0:7), nBas_Valence(0:7)
      Integer indexLarge,indexSmall,iIrrep

      If(nIrrep.eq.1) return

      indexLarge = nDens + 1
      indexSmall = nDens_Valence + 1
      Do iIrrep = nIrrep-1, 0, -1
! calculate the position in the hypothetical Array(nDens_Valence)
! and in the needed Array(nDens)
        indexLarge = indexLarge - nBas(iIrrep)
        indexSmall = indexSmall - nBas_Valence(iIrrep)
! move the data
        call dcopy_(nBas_Valence(iIrrep), Array(indexSmall), 1,         &
     &                                   Array(indexLarge), 1)
        call dcopy_(nBas_Valence(iIrrep),[Zero],0,Array(indexSmall),1)
      End Do
      Return
      End
