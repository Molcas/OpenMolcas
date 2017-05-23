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
      Subroutine LDF_SetThrs(Target_Accuracy)
      Implicit None
      Real*8 Target_Accuracy
#include "localdf.fh"
      Real*8 CutInt
      Thr_Accuracy=Target_Accuracy
      If (Thr_Prescreen.lt.0.0d0) Then
         Call LDF_GetCutInt(CutInt)
         Call LDF_SetPrescreen(min(Thr_Accuracy,CutInt))
      End If
      End
