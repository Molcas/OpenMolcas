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
      SUBROUTINE CHO_PTRINI(irc)
C
C     Purpose: set all entries in /CHMIND/ zero.
C
      IMPLICIT NONE
      Integer irc
#include "choptr.fh"
      Integer nAlloc

      nAlloc = 0  ! allocation counter

      ip_INFVEC = 0
      l_INFVEC  = 0
      nAlloc    = nAlloc + 1

      ip_INDRED = 0
      l_INDRED  = 0
      nAlloc    = nAlloc + 1

      irc = CHO_NALLOC - nAlloc

      End
