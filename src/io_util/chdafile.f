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
      Subroutine ChDaFile(Lu,iOpt,Buf,lBuf,iDisk)
      Implicit Integer (A-Z)
      Character Buf(*)

      Call ChDaFile_Internal(Buf)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine ChDaFile_Internal(Buf)
      Use Iso_C_Binding
      Character, Target :: Buf(*)
      Integer, Pointer :: iBuf(:)
      Call C_F_Pointer(C_Loc(Buf(1)),iBuf,[1])
      Call DaFile(Lu,iOpt,iBuf,lBuf,iDisk)
      Nullify(iBuf)
      End Subroutine ChDaFile_Internal
*
      End
