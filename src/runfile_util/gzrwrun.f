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
* Copyright (C) 2009, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine reads or writes a record from/to the runfile.           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: December 2009                                               *
*                                                                      *
************************************************************************
      Subroutine gzRWRun(Lu,icXX,Data,nData,iDisk,RecTyp)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       Lu
      Integer       icXX
      Character     Data(*)
      Integer       nData
      Integer       iDisk
      Integer       RecTyp
*----------------------------------------------------------------------*
* Read/write data from/to runfile.                                     *
*----------------------------------------------------------------------*
      If(RecTyp.eq.TypInt) Then
         Call c_iDaFile(Lu,icXX,Data,nData,iDisk)
      Else If(RecTyp.eq.TypDbl) Then
         Call c_dDaFile(Lu,icXX,Data,nData,iDisk)
      Else If(RecTyp.eq.TypStr) Then
         Call cDaFile(Lu,icXX,Data,nData,iDisk)
      Else If(RecTyp.eq.TypLgl) Then
         Call SysAbendMsg('gzRWRun',
     &                    'Records of logical type not implemented',
     &                    'Aborting')
      Else
         Call SysAbendMsg('gzRWRun',
     &                    'Argument RecTyp is of wrong type',
     &                    'Aborting')
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine c_iDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)
      Use Iso_C_Binding
      Integer Lu, iOpt, lBuf_, iDisk_
      Character, Target :: Buf(*)
      Integer, Pointer :: pBuf(:)
      Call C_F_Pointer(C_Loc(Buf(1)),pBuf,[lBuf_])
      Call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
      Nullify(pBuf)
      End SubRoutine c_iDaFile
      SubRoutine c_dDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)
      Use Iso_C_Binding
      Integer Lu, iOpt, lBuf_, iDisk_
      Character, Target :: Buf(*)
      Real*8, Pointer :: pBuf(:)
      Call C_F_Pointer(C_Loc(Buf(1)),pBuf,[lBuf_])
      Call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
      Nullify(pBuf)
      End SubRoutine c_dDaFile
*
      End
