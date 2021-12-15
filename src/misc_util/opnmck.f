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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine OpnMCK (rc,Option,Name,Lu)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Open the one-electron integral file and initialize/load          *
*     the table of contents.                                           *
*                                                                      *
*     input:                                                           *
*     option : Switch to set options                                   *
*              = 0 old file                                            *
*              = 1 new file                                            *
*     FnCom  : Logical file name                                       *
*     LuCom  : FORTRAN unit number                                     *
*                                                                      *
*     output:                                                          *
*     rc     : Return code.                                            *
*              A value of 0 (zero) is returned upon successful         *
*              completion of the request. A nonzero value indi-        *
*              cates an error.                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Integer (A-Z)
*
#include "FileIDs.fh"
#include "MckRc.fh"
#include "MckFlags.fh"
#include "SysDef.fh"
#include "MckDat.fh"

*
      Character*(*) Name
      Character*8   FnMCK
      Logical exist,NewToc
      Character*16 TheName
      Data TheName/'OpnMck'/
*---------------------------------------------------------------------*
*     Start procedure:                                                *
*---------------------------------------------------------------------*
      NewToc=iAnd(option,sNew).ne.0
      rc=rc0000
*---------------------------------------------------------------------*
*     Start procedure:                                                *
*---------------------------------------------------------------------*
      AuxMCK(pLu   ) = 0
      AuxMCK(pOpen ) = 0
      Call StdFmt(Name,FnMCK)
      LuMCK=Lu
      call f_Inquire ( FnMCK,Exist)
*----------------------------------------------------------------------*
*     Check the options                                                *
*----------------------------------------------------------------------*
      If ( Option.ne.0 ) Then
         SumOpt=0
         If ( iAnd(Option,sNew).ne.0 ) SumOpt=SumOpt+sNew
         If ( iAnd(Option,1024).ne.0 ) SumOpt=SumOpt+1024
         If ( SumOpt.ne.Option ) Then
           Call SysWarnMsg(TheName,'MSG: invalid option',' ')
       Call SysCondMsg('SumOpt.eq.Option',SumOpt,'<>',Option)
         End If
      End If
*----------------------------------------------------------------------*
*     Compare file status with options                                 *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Old file did not exist                                           *
*----------------------------------------------------------------------*
      If ( .not.NewToc .and. .not.exist) then
           Call SysAbendMsg(TheName,
     * 'MCK file does not exist',' ')
*----------------------------------------------------------------------*
*     New toc                                                          *
*----------------------------------------------------------------------*
      Else If( NewToc ) Then
         Call DaName(LuMCK,FnMCK)
         Call iCopy(lToc,[NaN],0,TocOne,1)
         TocOne(pFID)=IDone
         TocOne(pVersN)=VNone
         iDisk=0
         Call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
         TocOne(pNext)=iDisk
         iDisk=0
         Call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
         AuxMCK(pLu   ) = LuMCK
         AuxMCK(pOpen ) = 1
*----------------------------------------------------------------------*
*     Old toc                                                          *
*----------------------------------------------------------------------*
      Else
         Call DaName(LuMCK,FnMCK)
         iDisk=0
         Call iDaFile(LuMCK,2,TocOne,lToc,iDisk)
         If( TocOne(pFID).ne.IDone .or. TocOne(pVersN).ne.VNone ) Then
           Call SysFileMsg(TheName,
     * 'file version number is outdated',LuMCK,' ')
         End If
         AuxMCK(pLu   ) = LuMCK
         AuxMCK(pOpen ) = 1
      End If
*     Report back the actual unit number
      Lu=LuMCK
*----------------------------------------------------------------------*
*     normal end                                                       *
*----------------------------------------------------------------------*
      Return
      End
