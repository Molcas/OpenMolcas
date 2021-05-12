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
      Subroutine LDF_GetThrs()
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Get LDF target accuracy from Runfile
C
      Implicit None
      Real*8 Target_Accuracy
      Call Get_dScalar('LDF Accuracy',Target_Accuracy)
      Call LDF_SetThrs(Target_Accuracy)
#if defined (_DEBUGPRINT_)
      Write(6,'(A,1P,D15.6)')
     & 'LDF_GetThrs: target accuracy from runfile:',Target_Accuracy
#endif
      End
