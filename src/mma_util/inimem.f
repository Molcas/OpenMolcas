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
*  IniMem
*
*> @brief
*>   Initialize memory for Molcas
*>
*> @details
*> Initialize memory for Molcas.
************************************************************************
      Subroutine IniMem
      Use stdalloc, only: MxMem
      Implicit Real*8 (A-H,O-Z)
*
#include "SysCtl.fh"
#include "warnings.h"
#include "mama.fh"
#include "WrkSpc.fh"
*
      Interface
        Function allocmem(ref,intof,dblof,chrof,size_)
     &           bind(C,name='allocmem_')
          Use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
          Integer(kind=MOLCAS_C_INT) :: allocmem
          Real(kind=MOLCAS_C_REAL) :: ref(*)
          Integer(kind=MOLCAS_C_INT) :: intof, dblof, chrof, size_
        End Function allocmem
      End Interface
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
      MemCtl(ipSysOut) = 6

*----------------------------------------------------------------------*
*     Grab from the system a pointer to the dynamic work area          *
*----------------------------------------------------------------------*
      iRc=allocmem(Work,iofint,iofdbl,iofchr,MxMem)
      If ( iRc.ne.0 ) Then
         Write (6,'(A,I3,A)') 'The initialization of the memory '//
     &                        'manager failed ( iRc=',iRc,' ).'
         Call Quit(_RC_MEMORY_ERROR_)
      End If
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      End Subroutine
