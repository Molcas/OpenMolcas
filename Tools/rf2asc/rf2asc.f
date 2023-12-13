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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This program converts a runfile into an ascii file.                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
* Written: July 2003                                                   *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Program RF2Asc
      Use RunFile_data, Only: lw, nToc, NulPtr, RunHdr, Toc, TypDbl,    &
     &                        TypInt, TypStr
*----------------------------------------------------------------------*
* Declare local variables.                                             *
*----------------------------------------------------------------------*
      Integer   MxBuf
      Parameter ( MxBuf = 1024*1024 )
      Integer   Lu,iDisk
      Integer   i,j
      Real*8, Allocatable, Dimension(:)    :: dBuf
      Integer, Allocatable, Dimension(:)   :: iBuf
      Character, Allocatable, Dimension(:) :: cBuf
      Integer   iRc
      Integer   iOpt
      Call Init_LinAlg
      Call PrgmInit('RF2Asc')
*----------------------------------------------------------------------*
* Open runfile and check that file is ok.                              *
*----------------------------------------------------------------------*
      Call NameRun('RUNFILE')
      iOpt=0
      iRc=0
      Call OpnRun(iRc,Lu, iOpt)
*----------------------------------------------------------------------*
* Read the ToC                                                         *
*----------------------------------------------------------------------*
      iDisk=RunHdr%DaLab
      Call cDaFile(Lu,icRd,Toc(:)%Lab,lw*nToc,iDisk)
      iDisk=RunHdr%DaPtr
      Call iDaFile(Lu,icRd,Toc(:)%Ptr,nToc,iDisk)
      iDisk=RunHdr%DaLen
      Call iDaFile(Lu,icRd,Toc(:)%Len,nToc,iDisk)
      iDisk=RunHdr%DaTyp
      Call iDaFile(Lu,icRd,Toc(:)%Typ,nToc,iDisk)
*----------------------------------------------------------------------*
* Open output file.                                                    *
*----------------------------------------------------------------------*
      Open(Unit=9, File='RUNASCII', Err=999)
*----------------------------------------------------------------------*
* Print header information.                                            *
*----------------------------------------------------------------------*
      Write(9,'(a)')     '# Header information'
      Write(9,'(a,i15)') '# ID                      ',RunHdr%ID
      Write(9,'(a,i15)') '# Version                 ',RunHdr%Ver
      Write(9,'(a,i15)') '# Next free address       ',RunHdr%Next
      Write(9,'(a,i15)') '# Number of items         ',RunHdr%Items
      Write(9,'(a,i15)') '# Address to ToC labels   ',RunHdr%DaLab
      Write(9,'(a,i15)') '# Address to ToC pointers ',RunHdr%DaPtr
      Write(9,'(a,i15)') '# Address to ToC lengths  ',RunHdr%DaLen
      Write(9,'(a,i15)') '# Address to ToC types    ',RunHdr%DaTyp
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Allocate(dBuf(MxBuf),iBuf(MxBuf),cBuf(MxBuf))
      Do i=1,nToc
         If(Toc(i)%Ptr.ne.NulPtr) Then
            Write(9,'(3a)') '<',Toc(i)%Lab,'>'
            Write(9,*) Toc(i)%Len,Toc(i)%Typ
            If(Toc(i)%Typ.eq.TypDbl) Then
               iDisk=Toc(i)%Ptr
               Call dDaFile(Lu,icRd,dBuf,Toc(i)%Len,iDisk)
               Write(9,'(4d26.18)') (dBuf(j),j=1,Toc(i)%Len)
            Else If(Toc(i)%Typ.eq.TypInt) Then
               iDisk=Toc(i)%Ptr
               Call iDaFile(Lu,icRd,iBuf,Toc(i)%Len,iDisk)
               Write(9,*) (iBuf(j),j=1,Toc(i)%Len)
            Else If(Toc(i)%Typ.eq.TypStr) Then
               iDisk=Toc(i)%Ptr
               Call cDaFile(Lu,icRd,cBuf,Toc(i)%Len,iDisk)
               Write(9,'(64a1)') (cBuf(j),j=1,Toc(i)%Len)
            Else
            End If
         End If
      End Do
      Deallocate(dBuf,iBuf,cBuf)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Close(9)
      Stop
999   Continue
      Write(*,*) 'An error occured'
      Stop
      End
