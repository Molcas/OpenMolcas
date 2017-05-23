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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine ClsFls_CASPT2
************************************************************************
*     Close files.                                                     *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)
*----------------------------------------------------------------------*
*     Start                                                            *
*-------------------------------------- -------------------------------*

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
      Call qEnter('ClsFls')
      Call DaClos(LUCIEX)
* PAM08
*      Call DaClos(LUMORB)
      Call DaClos(LUONEM)
      Call DaClos(LUINTM)
      Call DaClos(LUDRA)
      Call DaClos(LUDRATOT)
      Call DaClos(LUHLF1)
      Call DaClos(LUHLF2)
      Call DaClos(LUHLF3)
      Call DaClos(LUDMAT)
      Call DaClos(LUSOLV)
      Call DaClos(LUSBT)
      DO IVEC=1,6
        CALL DaClos(LURHS(IVEC))
      END DO
      DO IMAT=1,4
        CALL DaClos(LUH0T(IMAT))
      END DO
*---  close the ORDINT file -------------------------------------------*
      If (.not.IfChol) then
         iRc=-1
         iOpt=0
         Call ClsOrd(iRc,iOpt)
         IF(IRC.NE.0 .AND. IPRGLB.GT.SILENT) THEN
          Call WarningMessage(1,'Failed to close ORDINT file.')
         END IF
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Call qExit('ClsFls')
      Return
      End
