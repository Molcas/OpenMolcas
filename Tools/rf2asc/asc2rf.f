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
      Program Asc2RF
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare local variables.                                             *
*----------------------------------------------------------------------*
      Integer      MxBuf
      Parameter    ( MxBuf = 1024*1024 )
      Real*8       dBuf(MxBuf)
      Integer      iBuf(MxBuf)
      Character    cBuf(MxBuf)

      Character*18 Rec
      Character*16 Label
      Integer      Length
      Integer      Type
      Integer      iRc
      Integer      iOpt
      Integer      i
      Call Init_LinAlg
*----------------------------------------------------------------------*
* Create RunFile.                                                      *
*----------------------------------------------------------------------*
      Call NameRun('RUNFILE')
      iOpt=0
      iRc=0
      Call MkRun(iRc, iOpt)
*----------------------------------------------------------------------*
* Open input file.                                                     *
*----------------------------------------------------------------------*
      Open(Unit=9, File='RUNASCII', Status='OLD', Err=999)
*----------------------------------------------------------------------*
* Populate runfile.                                                    *
*----------------------------------------------------------------------*
  100 Continue
      Read(9,'(a18)',End=998,Err=999) Rec
         if(Rec(1:1).eq.'#') GoTo 100
         If(Rec(1:1).ne.'<') GoTo 999
         If(Rec(18:18).ne.'>') GoTo 999
         Label=Rec(2:17)
         Read(9,*) Length,Type
         Write(*,*) 'Processing:',Rec,Length,Type


         If(Type.eq.TypInt) Then
            Read(9,*) (iBuf(i),i=1,Length)
            Call ixWrRun(iRc,Label, iBuf,Length, iOpt)
         Else If(Type.eq.TypDbl) Then
            Read(9,*) (dBuf(i),i=1,Length)
            Call dxWrRun(iRc,Label, dBuf,Length, iOpt)
         Else If(Type.eq.TypStr) Then
            Read(9,'(64a1)') (cBuf(i),i=1,Length)
            Call cxWrRun(iRc,Label, cBuf,Length, iOpt)
         Else If(Type.eq.TypLgl) Then
            Write(*,*) 'Cannot handle type logical'
            Stop
         Else
            Write(*,*) 'Unknown data type:',Type
         End If
         GoTo 100
      Continue
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
  998 Continue
      Close(9)
      Stop

  999 Continue
      Write(*,*) 'An error occured'
      Stop
      End
