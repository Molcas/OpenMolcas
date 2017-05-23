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
      Subroutine Put_SCF_Info_I(iCase,Arr,nArr)
      Implicit Real*8 (A-H,O-Z)

      Integer      Arr(nArr)
      Character*24 Label

      Label='SCFInfoI'
      if(iCase.eq.1) Label='SCFInfoI_ab'
      Call Put_iArray(Label,Arr,nArr)

      Return
      End
