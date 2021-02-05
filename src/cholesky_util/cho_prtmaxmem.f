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
      SubRoutine Cho_PrtMaxMem(Location)
C
C     Purpose: print max. available memory block.
C
      Implicit None
      Character(LEN=*) Location
#include "cholesky.fh"

      Character(LEN=2) Unt
      Integer l, lMax
      Real*8 dlMax

      l = len(Location)
      If (l .lt. 1) Then
         Write(Lupri,'(/,A)')
     &   'Largest memory block available @<UNKNOWN>:'
      Else
         Write(Lupri,'(/,A,A,A)')
     &   'Largest memory block available @',Location(1:l),':'
      End If
      Call mma_maxDBLE(lMax)
      Call Cho_Word2Byte(lMax,8,dlMax,Unt)
      Write(Lupri,'(3X,I10,A,F10.3,A,A)')
     & lMax,' 8-byte words; ',dlMax,' ',Unt

      End
