************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
************************************************************
*
*   <DOC>
*     <Name>IniMem</Name>
*     <Syntax>Call IniMem()</Syntax>
*     <Arguments>
*     </Arguments>
*     <Purpose>
* To initialize memory for the Molcas.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
* To initialize memory for the Molcas.
*     </Description>
*    </DOC>
*
************************************************************
      Subroutine IniMem
      Implicit Real*8 (A-H,O-Z)
*
#include "SysCtl.fh"
#include "warnings.fh"
#include "mama.fh"
#include "WrkSpc.fh"
*
      Integer allocmem
*----------------------------------------------------------------------*
*     Initialize the Common / MemCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      Do i=1,ipCheck
         MemCtl(i)=0
      End Do
      MemCtl(ipStat)   = ON
      MemCtl(ipTrace)  = OFF
      MemCtl(ipQuery)  = OFF
      MemCtl(ipCheck)  = OFF
      MemCtl(ipClear)  = OFF
#ifndef NAGFOR
#ifdef _GARBLE_
      MemCtl(ipCheck)  = ON
      MemCtl(ipClear)  = ON
#endif
#endif
      MemCtl(ipSysOut) = 6

*----------------------------------------------------------------------*
*     Grab from the system a pointer to the dynamic work area          *
*----------------------------------------------------------------------*
      iRc=allocmem(Work,cWork,iofint,iofdbl,iofsgl,iofchr,MxMem)
      If ( iRc.ne.0 ) Then
         Write (6,'(A,I3,A)') 'The initialization of the memory '//
     &                        'manager failed ( iRc=',iRc,' ).'
         Call Quit(_RC_MEMORY_ERROR_)
      End If
*----------------------------------------------------------------------*
*     Allocate "dummy" pointers                                        *
*----------------------------------------------------------------------*
      Call GetMem('ip_Dum', 'Allo','REAL',ip_Dummy,  1)
      Call GetMem('ip_sDum','Allo','SNGL',ip_sDummy, 1)
      Call GetMem('ip_iDum','Allo','INTE',ip_iDummy, 1)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
