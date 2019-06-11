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
      Subroutine dRdMCK(rc,Option,InLab,iComp,dData,iSymLab)
      Implicit Integer (A-Z)
      Character*(*) InLab
      Real*8 dData(*)
      Call dRdMCK_Internal(dData)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine dRdMCK_Internal(dData)
      Use Iso_C_Binding
      Real*8, Target :: dData(*)
      Integer, Pointer :: iData(:)
      Call C_F_Pointer(C_Loc(dData),iData,[1])
      Call RdMCK(rc,Option,InLab,iComp,iData,iSymLab)
      Nullify(iData)
      Return
      End Subroutine dRdMCK_Internal
*
      End
