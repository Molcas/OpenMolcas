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
************************************************************************
      Subroutine OpnOrd (rc,Option,Name,Lu)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Open the two-electron integral file and initialize/load          *
*     the table of contents.                                           *
*                                                                      *
*     input:                                                           *
*     option : Switch to set options                                   *
*              = 0 old file                                            *
*              = 1 new file                                            *
*     Name   : Logical file name                                       *
*     Lu     : FORTRAN unit number                                     *
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
*     M. P. Fuelscher                                                  *
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
#include "TwoRc.fh"
#include "TwoFlags.fh"
#include "Molcas.fh"
#include "TwoDat.fh"
*
      Parameter (NaN = -1)
      Character*(*) Name
      Character*8   FnTwo
      Logical Exist,NewToc
*
      Logical lDummy
      Integer nDummy1(8), nDummy2(8)
      Character*16 TheName
      Data TheName/'OpnOrd'/
**---------------------------------------------------------------------*
*     Start procedure:                                                *
*---------------------------------------------------------------------*
      rc=rc0000
*---------------------------------------------------------------------*
*     Check the file status                                           *
*---------------------------------------------------------------------*
      AuxTwo(isUnit) = iNoNum
      AuxTwo(isStat) = iNoNum
      AuxTwo(isDaDa) = iNoNum
      TocTwo(isPkPa) = iNoNum
      TocTwo(isPkAs) = iNoNum
      Call StdFmt(Name,FnTwo)
      LuTwo=Lu
      call f_Inquire ( FnTwo, Exist)
*----------------------------------------------------------------------*
*     Check the options                                                *
*----------------------------------------------------------------------*
      If ( Option.ne.0 ) Then
         SumOpt=0
         If ( iAnd(Option,sNew).ne.0 ) SumOpt=SumOpt+sNew
         If ( SumOpt.ne.Option ) Then
            Call SysWarnMsg(TheName,'MSG: invalid option',' ')
      Call SysCondMsg('SumOpt.eq.Option',SumOpt,'<>',Option)
         End If
      End If
      NewToc=iAnd(sNew,Option).ne.0
*----------------------------------------------------------------------*
*     Compare file status with options                                 *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Old file did not exist                                           *
*----------------------------------------------------------------------*
      If (.not.Exist.and..not.NewToc) Then
           Call SysAbendMsg(TheName,
     * 'ORDINT file does not exist',' ')
*----------------------------------------------------------------------*
*     New Toc                                                          *
*----------------------------------------------------------------------*
      Else If( NewToc )  Then
         Call DaName_MF(LuTwo,FnTwo)
         Call iCopy(lTocTwo,[NaN],0,TocTwo,1)
         TocTwo(isId)=IDtwo
         TocTwo(isVer)=VNtwo
         TocTwo(isForm)=0
         iDisk=0
         Call iDaFile(LuTwo,1,TocTwo,lTocTwo,iDisk)
         AuxTwo(isUnit) = LuTwo
         AuxTwo(isStat) = 1
         AuxTwo(isDaDa) = 0
*----------------------------------------------------------------------*
*     Keep Toc                                                         *
*----------------------------------------------------------------------*
      Else
         Call DaName_MF(LuTwo,FnTwo)
         iDisk=0
         Call iDaFile(LuTwo,2,TocTwo,lTocTwo,iDisk)
         If( TocTwo(isId).ne.IDtwo .or. TocTwo(isVer).ne.VNtwo ) Then
           Call SysFileMsg(TheName,
     * 'file version number is outdated',LuTwo,' ')
         End If
         AuxTwo(isUnit) = LuTwo
         AuxTwo(isStat) = 1
         AuxTwo(isDaDa) = iDisk
      End If
*
*---- Call to GetOrd to fill nBatch etc.
*
      If (Option.eq.0) Then
         rc_Dummy=-1
         Call GetOrd(rd_Dummy,lDummy,iDummy,nDummy1,nDummy2)
      End If
*
*----------------------------------------------------------------------*
*     normal end                                                       *
*----------------------------------------------------------------------*
      Return
      End
