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
C
C     Purpose: to deallocate memory of the  Cholesky MP2 program.
C
#include "implicit.fh"
#include "chomp2.fh"

      irc = 0

      Call ChoMP2g_deallocate(irc)

      If (.NOT.ChoMP2_allocated) Return

      Call GetMem('LiPQprod','Free','Inte',ip_LiPQprod,l_LiPQprod)
      Call GetMem('LnPQprod','Free','Inte',ip_LnPQprod,l_LnPQprod)
      Call GetMem('LnBatOrb','Free','Inte',ip_LnBatOrb,l_LnBatOrb)
      Call GetMem('NumBatOrb','Free','Inte',ip_NumBatOrb,l_NumBatOrb)
      Call GetMem('lUnit','Free','Inte',ip_lUnit,l_lUnit)
      Call GetMem('LiMatij','Free','Inte',ip_LiMatij,l_LiMatij)
      Call GetMem('LnMatij','Free','Inte',ip_LnMatij,l_LnMatij)
      Call GetMem('LiT1am','Free','Inte',ip_LiT1am,l_LiT1am)
      Call GetMem('LnT1am','Free','Inte',ip_LnT1am,l_LnT1am)
      Call GetMem('LnOcc','Free','Inte',ip_LnOcc,l_LnOcc)
      Call GetMem('NumOcc','Free','Inte',ip_NumOcc,l_NumOcc)
      Call GetMem('FirstS','Free','Inte',ip_FirstS,l_FirstS)
      Call GetMem('First','Free','Inte',ip_First,l_First)
*
      l_LiT1am  = 0
      l_LnT1am  = 0
      l_LnOcc   = 0
      l_NumOcc  = 0
      l_First   = 0
      l_FirstS  = 0
      l_LnMatij = 0
      l_LiMatij = 0
      l_lUnit   = 0
      l_NumBatOrb = 0
      l_LnBatOrb  = 0
      l_LnPQprod = 0
      l_LiPQprod = 0

      ChoMP2_allocated=.FALSE.

      End

      SubRoutine ChoMP2g_deallocate(irc)
      use ChoMP2, only: ChoMP2g_allocated
*
*     Purpose: Deallocate memory needed for
*              MP2-gradients or properties.
*
#include "implicit.fh"
#include "chomp2g.fh"
#include "chomp2.fh"

      irc = 0

      If (.NOT.ChoMP2g_allocated) Return

      Call GetMem('MoMoTable','Free','Inte',ipMoMoTable,lMoMoTable)
      Call GetMem('MP2Density','Free','Real',ipMP2D, lDens)
      Call GetMem('MP2WDensity','Free','Real',ipMP2W, lDens)
      Call GetMem('MP2Density_e','Free','Real',ipMP2D_e, lDens_e)
      Call GetMem('MP2WDensity_e','Free','Real',ipMP2W_e, lDens_e)
      Call GetMem('AdrVector1','Free','Inte',ipAdrR1, lAdrR1)
      Call GetMem('AdrVector2','Free','Inte',ipAdrR2, lAdrR2)
      Call GetMem('EFro','Free','Real',ip_EFroz,nFroT)
      Call GetMem('EOcc','Free','Real',ip_EOccu,nOccT)
      Call GetMem('EVir','Free','Real',ip_EVirt,nVirT)

      ChoMP2g_allocated=.FALSE.

      End
