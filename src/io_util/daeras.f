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
      Subroutine DaEras(Lu)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Close unit Lu and remove the file                                *
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
      Character*80 Text
      Character*16 TheName
      Data TheName/'DaEras'/

      If ( Query ) Call qEnter(TheName)

      If ( Trace ) then
        Write (6,*) ' >>> Enter DaEras <<<'
        Write (6,*) ' unit :',Lu
      End If

*     Save calling arguments

      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     * Call SysFileMsg(TheName,'MSG: unit', Lu,' ')

      If ( isOpen(Lu).eq.0 )
     * Call SysFileMsg(TheName,'MSG: used', Lu,' ')
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      If(isFiM(Lu).eq.0) then
#endif
        iRc = AixCls(FSCB(Lu))
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      Else
        iRc=FimCls(FSCB(Lu))
      End If
#endif
      If ( iRc.ne.0 ) then
        iRc = AixErr(Text)
        Call SysFileMsg(TheName,'MSG: close', Lu,Text)
      End If
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      If(isFiM(Lu).eq.0) then
#endif
        iRc = AixRm(LuName(Lu))
        If ( iRc.ne.0 ) then
          iRc = AixErr(Text)
        Call SysFileMsg(TheName,'MSG: delete', Lu,Text)
        End If
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      End If
#endif
      isOpen(Lu) = 0

      If ( Multi_File(Lu) ) then
        If ( MaxFileSize.ne.0 ) then
          If ( Trace ) Then
            Write (6,*) ' This is a partitioned data set'
          End If
          Do i = 1,MaxSplitFile-1
            Lu_=MPUnit(i,Lu)
            If ( Lu_.ge.1 ) Then
              If ( isOpen(Lu_).ne.0 ) Then
                iRc = AixCls(FSCB(Lu_))
                If ( iRc.ne.0 ) then
                  iRc = AixErr(Text)
                  Call SysFileMsg(TheName,'MSG: close', Lu_,Text)
                End If
                iRc = AixRm(LuName(Lu_))
                If ( iRc.ne.0 ) then
                  iRc = AixErr(Text)
                  Call SysFileMsg(TheName,'MSG: delete', Lu_,Text)
                End If
                isOpen(Lu_) = 0
              End If
            End If
          End Do
        End If
      End If

      If ( Trace ) then
        Write (6,*) ' >>> Exit DaEras <<<'
      End If

      If ( Query ) Call qExit(TheName)

      Return
      End
