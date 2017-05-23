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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************

      SubRoutine ChoMP2g_Density3(irc,CMO)
*     Jonas Bostrom, March 2010.
*
*     Purpose: Finalize MP2 Density.

      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(*)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
#include "chomp2g.fh"

      Character*8 ThisNm
      Character*16 SecNam
      Integer nOccAll(8), nOrbAll(8)
      Parameter (SecNam = 'ChoMP2g_Density3', ThisNm = 'Density3')

      Do iSym = 1, 8
         nOccAll(iSym) = nOcc(iSym) + nFro(iSym)
         nOrbAll(iSym) = nOrb(iSym) + nDel(iSym)
      End Do

      lTriDens = 0
      Do iSym = 1, nSym
         lTriDens = lTriDens + nOrbAll(iSym)*(nOrbAll(iSym)+1)/2
      End Do

      Do iSym = 1, nSym
         Do i = 1, nOrbAll(iSym)
            Do j = 1, nOrbAll(iSym)
               If((i.le. nOrb(iSym)) .and. (j.le. nOrb(iSym))) Then
                  Work(ipDensity_e(iSym) + i-1 + nOrbAll(iSym)*(j-1)) =
     &              Work(ipDensity(iSym) + i-1 + nOrb(iSym)*(j-1))
                  Work(ipWDensity_e(iSym)+ i-1 + nOrbAll(iSym)*(j-1)) =
     &             Work(ipWDensity(iSym) + i-1 + nOrb(iSym)*(j-1))
               Else
                  Work(ipDensity_e(iSym)+ i-1 + nOrbAll(iSym)*(j-1)) =
     &                     0.0d0
                  Work(ipWDensity_e(iSym)+ i-1 + nOrbAll(iSym)*(j-1)) =
     &                     0.0d0
               End If
            End Do
         End Do
      End Do
*
      Call GetMem('AOTriDens','Allo','Real',ipAOTriDens, lTriDens)
      Call GetMem('WAOTriDens','Allo','Real',ipWAOTriDens,lTriDens)
      Call FZero(Work(ipAOTriDens), lTriDens)
      Call FZero(Work(ipWAOTriDens), lTriDens)
*

      Call Build_Mp2Dens(ipAOTriDens,ipDensity_e,CMO,nSym,nOrbAll,
     &                   nOccAll,.True.)
      Call Build_Mp2Dens(ipWAOTriDens,ipWDensity_e,CMO,nSym,nOrbAll,
     &                   nOccAll,.False.)
      Call Put_D1ao_Var(Work(ipAOTriDens),lTriDens)
      Call Put_Fock_Occ(Work(ipWAOTriDens),lTriDens)
      Call GetMem('AOTriDens','Free','Real',
     &            ipAOTriDens, lTriDens)
      Call GetMem('WAOTriDens','Free','Real',
     &            ipWAOTriDens, lTriDens)
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(irc)
      End
