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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_Cleanup(irc)
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Clean up after RPA run (deallocate etc.)
C
      Implicit None
      Integer irc
#include "rpa_config.fh"
#include "rpa_data.fh"

      Integer  RPA_iUHF
      External RPA_iUHF

      Integer i

      irc=0

      ! Set "Relax Method" on Runfile
      Call Put_cArray('Relax Method',RPAModel,8)

      ! Deallocate memory
      Do i=1,RPA_iUHF()
         If (l_CMO(i).gt.0) Then
            Call GetMem('CMO(RPA)','Free','Real',ip_CMO(i),l_CMO(i))
         End If
         ip_CMO(i)=0
         l_CMO(i)=0
         If (l_EMO(i).gt.0) Then
            Call GetMem('EMO(RPA)','Free','Real',ip_EMO(i),l_EMO(i))
         End If
         ip_EMO(i)=0
         l_EMO(i)=0
         If (l_OccEn(i).gt.0) Then
            Call GetMem('OccEn','Free','Real',ip_OccEn(i),l_OccEn(i))
         End If
         ip_OccEn(i)=0
         l_OccEn(i)=0
         If (l_VirEn(i).gt.0) Then
            Call GetMem('OccEn','Free','Real',ip_VirEn(i),l_VirEn(i))
         End If
         ip_VirEn(i)=0
         l_VirEn(i)=0
      End Do

      End
