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
* Copyright (C) 1990, IBM                                              *
************************************************************************
      Integer Function iPrmt(jOper,iChct)
************************************************************************
*     Returns the phase factor of a basis function under a symmetry    *
*     operation, jOper. iChct contains the information about the       *
*     character of the basis function.                                 *
************************************************************************
      Use Symmetry_Info, only: iOper
      Implicit Real*8 (a-h,o-z)
#define _CHECK_
#ifdef _CHECK_
      If (Size(iOper)<1) Then
         Write (6,*) 'iPrmt; iOper not defined.'
         Call Abend()
      End If
#endif
      iPrmt = 1
      iCom= iAnd(iOper(jOper),iChct)
      Do i = 1, 3
         If (iAnd(iCom,2**(i-1)).ne.0) iPrmt = iPrmt*(-1)
      End Do
      Return
      End Function iPrmt
