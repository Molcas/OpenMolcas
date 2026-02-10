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
      Subroutine ClsFls_CASPT2()
************************************************************************
*     Close files.                                                     *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
************************************************************************
      use definitions, only: iwp
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: silent
      use caspt2_global, only: LUCIEX, LUONEM, LUHLF1, LUHLF2, LUHLF3,
     &                       LUINTM, LUDMAT, LUDRA, LUDRATOT, LURHS,
     &                       LUH0T, LUSOLV, LUSBT
      use caspt2_module, only: IfChol
      Implicit None
      integer(kind=iwp) IMAT, iRc, IVEC
*----------------------------------------------------------------------*
*     Start                                                            *
*-------------------------------------- -------------------------------*

      Call DaClos(LUCIEX)
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
      DO IVEC=1,8
        CALL DaClos(LURHS(IVEC))
      END DO
      DO IMAT=1,4
        CALL DaClos(LUH0T(IMAT))
      END DO
*---  close the ORDINT file -------------------------------------------*
      If (.not.IfChol) then
         iRc=-1
         Call ClsOrd(iRc)
         IF(IRC.NE.0 .AND. IPRGLB.GT.SILENT) THEN
          Call WarningMessage(1,'Failed to close ORDINT file.')
         END IF
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      End Subroutine ClsFls_CASPT2
