************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************

#ifndef _HAVE_EXTRA_

* Broadcast a file from the master to the slaves

      Subroutine PFGet_ASCII(FName)
#ifdef _MOLCAS_MPP_
      Use MPI
#endif
      Implicit None
      Character (Len=*), Intent(In) :: FName
#ifdef _MOLCAS_MPP_
#include "para_info.fh"
#include "mpp_info.fh"
#ifdef _I8_
      Integer*4, Parameter :: iType=MPI_INTEGER8
#else
      Integer*4, Parameter :: iType=MPI_INTEGER4
#endif
      Integer, Parameter :: LBuf=4096
      Character (Len=LBuf) :: Buf
      Integer :: LU, Err, FLen, Pos, Num
      Logical :: Failed

      ! Open the file for reading or writing
      ! Note that each process opens only one file, so there is a single
      ! unit number LU
      If (King()) Then
        Call Molcas_Open_Ext2(LU, FName, "stream", "unformatted", Err,
     &                        .False., 0, "old", Failed)
        If (Failed .or. (Err .ne. 0)) Then
          Write(6,*) "Failed to open file ", Trim(FName)
          Call AbEnd()
        End If
        Inquire(LU, Size=FLen)
      Else
        Call Molcas_Open_Ext2(LU, FName, "stream", "unformatted", Err,
     &                        .False., 0, "replace", Failed)
        If (Failed .or. (Err .ne. 0)) Then
          Write(6,*) "Failed to open file ", Trim(FName)
          Call AbEnd()
        End If
      End If
      ! Broadcast the file size
      Call MPI_BCAST(FLen, 1, iType, 0, MPI_COMM_WORLD, Err)
      If (Err .ne. 0) Then
        Write(6,*) "Failed to broadcast file size: ", FLen
        Call AbEnd()
      End If
      ! Pass the file content in chunks
      Pos=0
      Do While (Pos .lt. FLen)
        ! Length of this chunk
        Num = Min(LBuf, FLen-Pos)
        ! The master reads the file
        If (King()) Then
          Read(LU, IOStat=Err) Buf(1:Num)
          If (Err .ne. 0) Then
            Write(6,*) "Error reading the file ", Trim(FName)
            Call AbEnd()
          End If
        End If
        Call MPI_BCast(Buf(1:Num), Num, MPI_CHARACTER, mpp_rootid,
     &                 MPI_COMM_WORLD, Err)
        If (Err .ne. 0) Then
          Write(6,*) "Failed to broadcast message of length: ", Num
          Call AbEnd()
        End If
        ! The slaves write the file
        If (.not. King()) Then
          Write(LU, IOStat=Err) Buf(1:Num)
          If (Err .ne. 0) Then
            Write(6,*) "Error writing the file ", Trim(FName)
            Call AbEnd()
          End If
        End If
        Pos = Pos + LBuf
      End Do
      Close(LU)
      Return
#else
      ! Avoid unused argument warnings
      If (.False.) Call Unused_Character(FName)
#endif

      End Subroutine PFGet_ASCII

#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine Empty_PFGet_ASCII
      End Subroutine Empty_PFGet_ASCII
#endif
