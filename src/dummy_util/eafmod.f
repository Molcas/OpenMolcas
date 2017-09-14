#ifndef _HAVE_EXTRA_

      Subroutine EAFOpen(Lu, FName)
      Integer :: Lu
      Character(Len=*) :: FName
      Call Unused_Integer(Lu)
      Call Unused_Character(FName)
      End Subroutine EAFOpen

      Subroutine EAFClose(Lu)
      Integer :: Lu
      Call Unused_Integer(Lu)
      End Subroutine EAFClose

      Subroutine EAFAWrite(Lu, Buf, nBuf, Disk, id)
      Integer :: Lu, nBuf, Buf(nBuf), id
      Real*8 :: Disk
      Call Unused_Integer(Lu)
      Call Unused_Integer_Array(Buf)
      Call Unused_Integer(nBuf)
      Call Unused_Real(Disk)
      Call Unused_Integer(id)
      End Subroutine EAFAWrite

      Subroutine EAFARead(Lu, Buf, nBuf, Disk, id)
      Integer :: Lu, nBuf, Buf(nBuf), id
      Real*8 :: Disk
      Call Unused_Integer(Lu)
      Call Unused_Integer_Array(Buf)
      Call Unused_Integer(nBuf)
      Call Unused_Real(Disk)
      Call Unused_Integer(id)
      End Subroutine EAFARead

      Subroutine EAFWrite(Lu, Buf, nBuf, Disk)
      Integer :: Lu, nBuf, Buf(nBuf)
      Real*8 :: Disk
      Call Unused_Integer(Lu)
      Call Unused_Integer_Array(Buf)
      Call Unused_Integer(nBuf)
      Call Unused_Real(Disk)
      End Subroutine EAFWrite

      Subroutine EAFRead(Lu, Buf, nBuf, Disk)
      Integer :: Lu, nBuf, Buf(nBuf)
      Real*8 :: Disk
      Call Unused_Integer(Lu)
      Call Unused_Integer_Array(Buf)
      Call Unused_Integer(nBuf)
      Call Unused_Real(Disk)
      End Subroutine EAFRead

      Subroutine EAFWait(LU, id)
      Integer :: LU, id
      Call Unused_Integer(Lu)
      Call Unused_Integer(id)
      End Subroutine EAFWait

      Function Disk2Byte(Disk)
      Real*8 Disk2Byte, Disk
      Call Unused_Real(Disk)
      Call Unused_Real(Disk2Byte)
      End Function Disk2Byte

#endif
