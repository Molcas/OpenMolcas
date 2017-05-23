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
      SubRoutine Cho_PTS_Final(NVT,l_NVT)
C
C     Thomas Bondo Pedersen, April 2010.
C
      Implicit None
      Integer l_NVT
      Integer NVT(l_NVT)
#include "choglob.fh"

      Integer i

      Call iCopy(l_NVT,NVT,1,NumCho_G,1)
      NumChT_G=NumCho_G(1)
      Do i=2,l_NVT
         NumChT_G=NumChT_G+NumCho_G(i)
      End Do
      Call Cho_Final(.False.)

      End
