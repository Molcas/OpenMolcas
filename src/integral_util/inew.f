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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Integer Function iNew(iTest,iIrrep)
      use Symmetry_Info
      Integer iTest(8)
      iNew = 0
*-----Test iTest against all rows thus far.
      Do 10 i = 1, iIrrep
*--------Do a intger inner product.
         iGo = 0
         Do 11 j = 1, nIrrep
            iGo = iGo + iTest(j)*iChTbl(i-1,j-1)
 11      Continue
         If (iGo.ne.0) Then
*-----------Here if row already defined.
            iNew = i
            Return
         End If
 10   Continue
      iNew = iIrrep+1
      Return
      End
