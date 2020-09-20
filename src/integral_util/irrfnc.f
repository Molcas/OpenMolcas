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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      Integer Function IrrFnc(iFnc)
      use Symmetry_Info
      Integer iTest(8)
      ix = iAnd(iFnc,1)
      iy = iAnd(iFnc,2)/2
      iz = iAnd(iFnc,4)/4
      Do i = 1, nIrrep
         jx = iAnd(iOper(i-1),1)
         jy = iAnd(iOper(i-1),2)/2
         jz = iAnd(iOper(i-1),4)/4
         iCh = 1
         If (ix.ne.0 .and. jx.ne.0) iCh = -iCh
         If (iy.ne.0 .and. jy.ne.0) iCh = -iCh
         If (iz.ne.0 .and. jz.ne.0) iCh = -iCh
         iTest(i)=iCh
      End Do
      IrrFnc=iNew(iTest,nIrrep)-1
      Return
      End
