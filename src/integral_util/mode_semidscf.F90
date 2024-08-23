!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
!               1990, IBM                                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************
      Subroutine Mode_SemiDSCF(Wr_Mode)
      use IOBUF, only: iStatIO, Mode_Read, Disk, Disk_2, Mode_Write
      Implicit None
      Logical Wr_Mode
!
!     Write (6,*) 'Mode_SemiDSCF: Wr_Mode=',Wr_Mode
      If (Wr_Mode) Then
         If (iStatIO.eq.Mode_Read) Then
            Disk = Disk_2
            iStatIO = Mode_Write
!           Write (6,*) 'Changing to Write mode @',Disk
         End If
      Else
         If (iStatIO.eq.Mode_Write) Then
            Write (6,*) 'Change from Write to Read mode not implemented'
            Call Abend()
         End If
      End If
!
      End Subroutine Mode_SemiDSCF
