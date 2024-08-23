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
      Subroutine Init_SemiDSCF(FstItr,Thize,Cutint)
      use dEAF
      use IOBUF, only: IODone, Disk, iBuf, ipos, iStatIO, Mode_Write,   &
     &                 OnDisk, Mode_Read, Disk_1, Disk_2, lBuf, nBuf,   &
     &                 Buffer, ID, LuTmp
      implicit None
#include "SysDef.fh"
      Logical FstItr
      Real*8 Thize, CutInt

      Integer lBufOld, nBufOld
      Real*8 ThizeOld, CutIntOld
      real*8 control(4)
!     Write (6,*) 'Enter: Init_SemiDSCF'
!     Write (6,*) 'Ondisk=',Ondisk
!     Write (6,*) 'lBuf=',lBuf
!
!---- Initiate asynchronous double buffer I/O.
!
      IODone = .False.
      Disk = 0.0D0
      iBuf=1
      iPos = 1
      If (FstItr) Then
         iStatIO = Mode_Write
!        write(6,*) 'write istatio=',istatio
         control(1)=Dble(lbuf)
         control(2)=Dble(nbuf)
         control(3)=thize
         control(4)=cutint
!        write(6,*) 'control written:',control
!        Write (6,*) ' Initiate write @', Disk,'iBuf=',iBuf
         If(OnDisk) Call dEAFAwrite(LuTmp,control,4*RtoI,Disk,id)
      Else
         iStatIO = Mode_Read
!        write(6,*) 'read istatio=',istatio
!
!------- Initiate first read ahead of time.
!
!        Write (6,*) 'lBuf*RtoI=',lbuf*RtoI,' rtoi=',Rtoi
         If (OnDisk) then
!           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
            Call dEAFread(LuTmp,control,4*RtoI,Disk)
            Disk_2 = Disk
            Disk_1 = Disk
!           write(6,*) 'control read:',control
            lbufold=nint(control(1))
            nbufold=nint(control(2))
            thizeold=control(3)
            cutintold=control(4)
            if (lbufold.lt.lbuf) then
              write(6,*) 'Reducing the buffer size from ',lbuf,         &
     &                  ' to',lbufold
              lbuf=lbufold
            else if(lbufold.gt.lbuf) then
              write(6,*) 'Inconsistent buffer lengths. Old:',lbufold,   &
     &                   '  current:',lbuf
              call Abend()
            end if
            if(nbuf.ne.nbufold) then
              write(6,*) 'Inconsistent buffer number. Old:',nbufold,    &
     &                   '  current:',nbuf
              call Abend()
            end if
            if(abs(thize-thizeold).gt.1.d-10) then
              write(6,*) 'Resetting thize from',thize,' to',thizeold
              thize=thizeold
            end if
            if(cutintold.gt.cutint) then
              write(6,*) 'Inconsistent Cutint. Old:',cutintold,         &
     &                   '  current:',cutint
              call Abend()
            end if
!           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
!           If(OnDisk) Write (6,*) ' Initial EAFARead'
            Call dEAFARead(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk,id)
         End If
      End If
!
!     Write (*,*) 'Exit: Init_SemiDSCF'
      End Subroutine Init_SemiDSCF
