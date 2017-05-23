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
      Subroutine Qpg_zArray(Label,Found,nData)
      Implicit None

      Character*(*) Label
      Logical       FoundR, FoundI, Found
      Integer       nData, nDataR, nDataI

      call qpg_darray ('R'//Label,FoundR,nDataR)
      call qpg_darray ('I'//Label,FoundI,nDataI)

      if ((nDataR.eq.nDataI).and.FoundR.and.FoundI) then
         nData=nDataR
         found=.true.
      else
         found=.false.
         ndata=0
      end if

      return
      end
