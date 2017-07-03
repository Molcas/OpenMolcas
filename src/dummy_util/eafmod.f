#ifndef _HAVE_EXTRA_

      Subroutine EAFOpen(Lu, FName)
      Integer :: Lu
      Character(Len=*) :: FName
      End Subroutine EAFOpen

      Subroutine EAFClose(Lu)
      Integer :: Lu
      End Subroutine EAFClose

      Subroutine EAFAWrite(Lu, Buf, nBuf, Disk, id)
      Integer :: Lu, nBuf, Buf(nBuf), id
      Real*8 :: Disk
      End Subroutine EAFAWrite

      Subroutine EAFARead(Lu, Buf, nBuf, Disk, id)
      Integer :: Lu, nBuf, Buf(nBuf), id
      Real*8 :: Disk
      End Subroutine EAFARead

      Subroutine EAFWrite(Lu, Buf, nBuf, Disk)
      Integer :: Lu, nBuf, Buf(nBuf)
      Real*8 :: Disk
      End Subroutine EAFWrite

      Subroutine EAFRead(Lu, Buf, nBuf, Disk)
      Integer :: Lu, nBuf, Buf(nBuf)
      Real*8 :: Disk
      End Subroutine EAFRead

      Subroutine EAFWait(LU, id)
      Integer :: LU, id
      End Subroutine EAFWait

      Function Disk2Byte(Disk)
      Real*8 Disk2Byte, Disk
      End Function Disk2Byte

#endif
