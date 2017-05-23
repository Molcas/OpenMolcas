************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_P_DskSwp()
C
C     Purpose: swap pointers to local and global arrays InfRed and
C              InfVec.
C
C     NB: do NOT use along with Cho_P_IndxSwp, as the same swapping is
C         done there, too. Thus, calling both will have no net effect!
C
      Implicit None
#include "choptr.fh"
#include "choglob.fh"

      Integer iTmp

      iTmp = ip_InfRed_G
      ip_InfRed_G = ip_InfRed
      ip_InfRed = iTmp
      iTmp = l_InfRed_G
      l_InfRed_G = l_InfRed
      l_InfRed = iTmp

      iTmp = ip_InfVec_G
      ip_InfVec_G = ip_InfVec
      ip_InfVec = iTmp
      iTmp = l_InfVec_G
      l_InfVec_G = l_InfVec
      l_InfVec = iTmp

      End
