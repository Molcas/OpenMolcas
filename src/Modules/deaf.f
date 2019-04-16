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

* This module simply provides an interface for the EAF subroutines with
* real buffers instead of integers, note that nBuf is still the length
* of the integer buffer.

      Module dEAF
      Use Iso_C_Binding
      Implicit None

      Contains

      Subroutine dEAFARead(Lu,Buf,nBuf,Disk,id)
      Integer :: Lu, nBuf, id
      Real*8, Target :: Buf(nBuf)
      Real*8 :: Disk
      Integer, Pointer :: iBuf(:)
      Call C_F_Pointer(C_Loc(Buf),iBuf,[nBuf])
      Call EAFARead(Lu,iBuf,nBuf,Disk,id)
      Nullify(iBuf)
      End Subroutine dEAFARead

      Subroutine dEAFAWrite(Lu,Buf,nBuf,Disk,id)
      Integer :: Lu, nBuf, id
      Real*8, Target :: Buf(nBuf)
      Real*8 :: Disk
      Integer, Pointer :: iBuf(:)
      Call C_F_Pointer(C_Loc(Buf),iBuf,[nBuf])
      Call EAFAWrite(Lu,iBuf,nBuf,Disk,id)
      Nullify(iBuf)
      End Subroutine dEAFAWrite

      Subroutine dEAFRead(Lu,Buf,nBuf,Disk)
      Integer :: Lu, nBuf
      Real*8, Target :: Buf(nBuf)
      Real*8 :: Disk
      Integer, Pointer :: iBuf(:)
      Call C_F_Pointer(C_Loc(Buf),iBuf,[nBuf])
      Call EAFRead(Lu,iBuf,nBuf,Disk)
      Nullify(iBuf)
      End Subroutine dEAFRead

      Subroutine dEAFWrite(Lu,Buf,nBuf,Disk)
      Integer :: Lu, nBuf
      Real*8, Target :: Buf(nBuf)
      Real*8 :: Disk
      Integer, Pointer :: iBuf(:)
      Call C_F_Pointer(C_Loc(Buf),iBuf,[nBuf])
      Call EAFWrite(Lu,iBuf,nBuf,Disk)
      Nullify(iBuf)
      End Subroutine dEAFWrite

      End Module dEAF
