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
      SubRoutine Cho_X_GetIP_InfVec(InfVcT)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: get pointer to InfVec array for all vectors.
C
      use ChoSwp, only: InfVec, InfVec_Bak
      Implicit None
      Integer, Pointer:: InfVct(:,:,:)
#include "chopar.fh"
#if defined (_MOLCAS_MPP_)
#include "cho_para_info.fh"
#else
      Logical Cho_Real_Par
      Cho_Real_Par=.False.
#endif

      If (Cho_Real_Par) Then
         If (Allocated(InfVec_Bak)) Then
            InfVcT => InfVec_Bak
          Else
            Call Cho_Quit(
     &               'Initialization problem in Cho_X_GetIP_InfVec',103)
         End If
      Else
         InfVcT => InfVec
      End If

      End
