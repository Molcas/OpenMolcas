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
* Copyright (C) 1991, Per-Olof Widmark                                 *
*               1993,1996,1997, Markus P. Fuelscher                    *
*               1996, Luis Serrano-Andres                              *
*               2012, Victor P. Vysotskiy                              *
************************************************************************
      Subroutine DaClos(Lu)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Close unit Lu                                                    *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

#include "fio.fh"
#ifndef _OLD_IO_STAT_
#include "pfio.fh"
#endif

      Character*80 Text
      Character*16 TheName

      Data TheName/'DaClos'/

      If ( Query ) Call qEnter(TheName)

      If ( Trace ) then
        Write (6,*) ' >>> Enter DaClos <<<'
        Write (6,*) ' unit :',Lu
        Write (6,*) ' name :',LuName(Lu)
      End If
#ifndef _OLD_IO_STAT_
      LuP=0
      Do i=1,NProfFiles
         If(LuNameProf(i).eq.LuName(Lu)) Then
            LuP=i
         End If
      End Do
#ifdef _GA_
      FlsSize(LuP)=AixFsz(FSCB(Lu))
#else
      If(isFiM(Lu).eq.0) then
         FlsSize(LuP)=AixFsz(FSCB(Lu))
      Else
         FlsSize(LuP)=FimFsz(FSCB(Lu))
      End If
#endif
#endif

      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     * Call SysFileMsg(TheName,'MSG: unit', Lu,' ')
      If ( isOpen(Lu).eq.0 )
     * Call SysFileMsg(TheName,'MSG: notopened', Lu,' ')
#ifdef _GA_
       iRc = AixCls(FSCB(Lu))
#else
      If(isFiM(Lu).eq.0) then
       iRc = AixCls(FSCB(Lu))
      Else
        iRc=FimCls(FSCB(Lu))
        isFiM(Lu)=0
      End If
#endif
      If ( iRc.ne.0 ) then
        iRc = AixErr(Text)
      Call SysFileMsg(TheName,'MSG: close', Lu,Text)
      End If
      isOpen(Lu) = 0
      MBL(Lu)=0
      If ( Multi_File(Lu) ) then
        If ( MaxFileSize.ne.0 ) then
          If ( Trace ) Write (6,*) ' This is a partitioned data set'
          Do i = 1,MaxSplitFile-1
             Lu_=MPUnit(i,Lu)
             If (Lu_.ge.1) Then
                irc=0
                If ( isOpen(Lu_).ne.0 ) iRc = AixCls(FSCB(Lu_))
                If ( iRc.ne.0 ) then
                   iRc = AixErr(Text)
                   Call SysFileMsg(TheName,'MSG: close', Lu_,Text)
                End If
                isOpen(Lu_) = 0
                MPUnit(i,Lu)=-99
                Multi_File(Lu_)=.False.
                MBL(Lu_)=0
             End If
          End Do
        End If
        MPUnit(0,Lu)=0
        Multi_File(Lu)=.False.
      End If

      If ( Trace ) then
        Write (6,*) ' >>> Exit DaClos <<<'
      End If

      If ( Query ) Call qExit(TheName)

      Return
      End
