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
      Subroutine OpnOne (rc,Option,Name,Lu)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Open the one electron integral file                              *
*                                                                      *
*     input:                                                           *
*     option : Switch to set options                                   *
*     FnOne  : Logical file name                                       *
*     LuOne  : FORTRAN unit number                                     *
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
#include "OneRc.fh"
#include "OneFlags.fh"

#include "OneDat.fh"
*
      Character*(*) Name
      Character*8   FnOne
      Logical       Exist,NewToc
      Character*16 TheName
      Data TheName/'OpnOne'/
*---------------------------------------------------------------------*
*     Start procedure:                                                *
*---------------------------------------------------------------------*
*     If ( Query ) Call qEnter(TheName)
      rc=rc0000
*---------------------------------------------------------------------*
*     Get basis sets dimensions                                       *
*---------------------------------------------------------------------*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
*---------------------------------------------------------------------*
*     Truncate the name to 8 characters and convert it to upper case  *
*---------------------------------------------------------------------*
      LuOne=Lu
      FnOne=Name
      Call UpCase(FnOne)
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
      call f_Inquire ( FnOne,Exist)
      NewToc=iAnd(Option,sNew).ne.0
*----------------------------------------------------------------------*
*     Compare file status with options                                 *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Old file did not exist                                           *
      If ( .not.Exist.and..not.NewToc) then
           Call SysAbendMsg(TheName,
     * 'The ONEINT file does not exist',' ')
*----------------------------------------------------------------------*
*     New toc                                                          *
*----------------------------------------------------------------------*
       Else If (NewToc) Then
         Call iCopy(lAux,[NaN],0,AuxOne,1)
         Call iCopy(lToc,[NaN],0,TocOne,1)
         Call DaName_MF(LuOne,FnOne)
         TocOne(pFID)=IDrlx
         TocOne(pVersN)=VNrlx
         iDisk=0
         Call iDaFile(LuOne,sWrite,TocOne,lToc,iDisk)
         TocOne(pNext)=iDisk
         iDisk=0
         Call iDaFile(LuOne,sWrite,TocOne,lToc,iDisk)
         AuxOne(pLu   ) = LuOne
         AuxOne(pOpen ) = 1
       Else
*----------------------------------------------------------------------*
*     Keep toc                                                         *
*----------------------------------------------------------------------*
         Call DaName_MF(LuOne,FnOne)
         iDisk=0
         Call iDaFile(LuOne,sRead,TocOne,lToc,iDisk)
         If( TocOne(pFID).ne.IDrlx .or. TocOne(pVersN).ne.VNrlx ) then
           Call SysFileMsg(TheName,
     * 'file version number is outdated',LuOne,' ')
         End If
         AuxOne(pLu   ) = LuOne
         AuxOne(pOpen ) = 1
      End If
*----------------------------------------------------------------------*
*     Dump the TOC upon request                                        *
*----------------------------------------------------------------------*
      If ( iAnd(Option,1024).ne.0 ) Call DmpOne
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
*      If ( Query ) Call qExit(TheName)
      Return
      End
