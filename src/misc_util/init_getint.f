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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************

**************************************************************
* Initializes info needed by the Cholesky integral generator
*
* F. Aquilante
**************************************************************
      SUBROUTINE INIT_GETINT(RC)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER RC

#include "RdOrd.fh"

      rc = 0
      Call get_iscalar('nSym',nSym)
      Call get_iarray('nBas',nBas,nSym)
      Call INIT_NumCV(NumCho,nSym)

      LuCVec(1) = -1
      LuCVec(2) = -1

      pq1 = 0


      RETURN
      END

**************************************************************
**************************************************************
      SUBROUTINE INIT_NumCV(NumCV,nSymm)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER nSymm, NumCV(nSymm)

#include "cholesky.fh"

      Do i=1,nSymm
         NumCV(i)=NumCho(i)
      End Do

      RETURN
      END
