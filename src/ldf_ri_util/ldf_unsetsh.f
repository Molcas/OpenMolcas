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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_UnsetSh(irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Unset data in localdf_bas.fh
C
      Implicit None
      Integer irc
#include "localdf_bas.fh"

      ! set return code
      irc=0

      ! Deallocate iShlSO
      Call GetMem('LDF_iShlSO','Free','Inte',ip_iShlSO,l_iShlSO)

      ! Deallocate nBasSh
      Call GetMem('LDF_nBasSh','Free','Inte',ip_nBasSh,l_nBasSh)

      ! Deallocate iSOShl
      Call GetMem('LDF_iSOShl','Free','Inte',ip_iSOShl,l_iSOShl)

      ! Zero variables
      nBas_Auxiliary=0
      nBas_Valence=0
      nShell_Auxiliary=0
      nShell_Valence=0

      End
