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
C complex number in runfile
      Subroutine Get_zArray(Label,Data,nData)
      Implicit None

      Character*(*) Label
      Integer       nData
      real*8        RData(nData),IData(nData)
      Complex*16    Data(nData)

      call Get_dArray('R'//Label, RData, nData)
      call Get_dArray('I'//Label, IData, nData)

      Data=RData+(0.0,1.0)*IData

      Return
      End
