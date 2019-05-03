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
      Function C2R8(CBuf)
      Real*8 C2R8
      Character CBuf(*)
      C2R8=C2R8_Internal(CBuf)
*
*     This is to allow type punning without an explicit interface
      Contains
      Real*8 Function C2R8_Internal(CBuf)
      Use Iso_C_Binding
      Character, Target :: CBuf(*)
      Real*8, Pointer :: Buf
      Call C_F_Pointer(C_Loc(CBuf(1)),Buf)
      C2R8_Internal=Buf
      Nullify(Buf)
      Return
      End Function C2R8_Internal
*
      End
