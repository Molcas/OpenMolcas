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
      Integer Function iCLast(KWord,iChrct)
      Character KWord*(*)
      Do 10 i = iChrct, 1, -1
         If (KWord(i:i).ne.' ') Go To 11
 10   Continue
      iCLast = 0
      Return
 11   Continue
      iCLast = i
      Return
      End
