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
      Subroutine cRdMCK(rc,Option,InLab,iComp,cData,iSymLab)
      Implicit Integer (A-Z)
      Character*(*) InLab, cData
      Call cRdMCK_Internal(cData)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine cRdMCK_Internal(cData)
      Use Iso_C_Binding
      Character, Target :: cData(*)
      Integer, Pointer :: iData(:)
      Call C_F_Pointer(C_Loc(cData(1)),iData,[1])
      Call RdMCK(rc,Option,InLab,iComp,iData,iSymLab)
      Nullify(iData)
      Return
      End Subroutine cRdMCK_Internal
*
      End
