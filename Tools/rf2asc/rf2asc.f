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
#include "runinfo.fh"
#include "runtypes.fh"
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
      iDisk=RunHdr(ipDaLab)
      Call cDaFile(Lu,icRd,TocLab,16*nToc,iDisk)
      iDisk=RunHdr(ipDaPtr)
      Call iDaFile(Lu,icRd,TocPtr,nToc,iDisk)
      iDisk=RunHdr(ipDaLen)
      Call iDaFile(Lu,icRd,TocLen,nToc,iDisk)
      iDisk=RunHdr(ipDaTyp)
      Call iDaFile(Lu,icRd,TocTyp,nToc,iDisk)
*----------------------------------------------------------------------*
* Open output file.                                                    *
*----------------------------------------------------------------------*
      Open(Unit=9, File='RUNASCII', Err=999)
*----------------------------------------------------------------------*
* Print header information.                                            *
*----------------------------------------------------------------------*
      Write(9,'(a)')     '# Header information'
      Write(9,'(a,i15)') '# ID                      ',RunHdr(ipID)
      Write(9,'(a,i15)') '# Version                 ',RunHdr(ipVer)
      Write(9,'(a,i15)') '# Next free address       ',RunHdr(ipNext)
      Write(9,'(a,i15)') '# Number of items         ',RunHdr(ipItems)
      Write(9,'(a,i15)') '# Address to ToC labels   ',RunHdr(ipDaLab)
      Write(9,'(a,i15)') '# Address to ToC pointers ',RunHdr(ipDaPtr)
      Write(9,'(a,i15)') '# Address to ToC lengths  ',RunHdr(ipDaLen)
      Write(9,'(a,i15)') '# Address to ToC types    ',RunHdr(ipDaTyp)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Allocate(dBuf(MxBuf),iBuf(MxBuf),cBuf(MxBuf))
      Do i=1,nToc
         If(TocPtr(i).ne.NulPtr) Then
            Write(9,'(3a)') '<',TocLab(i),'>'
            Write(9,*) TocLen(i),TocTyp(i)
            If(TocTyp(i).eq.TypDbl) Then
               iDisk=TocPtr(i)
               Call dDaFile(Lu,icRd,dBuf,TocLen(i),iDisk)
               Write(9,'(4d26.18)') (dBuf(j),j=1,TocLen(i))
            Else If(TocTyp(i).eq.TypInt) Then
               iDisk=TocPtr(i)
               Call iDaFile(Lu,icRd,iBuf,TocLen(i),iDisk)
               Write(9,*) (iBuf(j),j=1,TocLen(i))
            Else If(TocTyp(i).eq.TypStr) Then
               iDisk=TocPtr(i)
               Call cDaFile(Lu,icRd,cBuf,TocLen(i),iDisk)
               Write(9,'(64a1)') (cBuf(j),j=1,TocLen(i))
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
