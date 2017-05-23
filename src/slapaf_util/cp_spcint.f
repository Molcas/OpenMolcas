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
      Subroutine cp_SpcInt
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Character*14  qLbl,filnam*16
      Logical Exist
*     Write (*,*) ' Copy files'
*
*.... Copy file SPCINX to SPCIN1
*
      filnam='SPCINX'
      Call f_Inquire(filnam,Exist)
      If (Exist) Then
      LuTmp1=11
      LuTmp2=12
      call molcas_binaryopen_vanilla(luTmp1,filnam)
c      Open(luTmp1,File=filnam,Form='unformatted',Status='unknown')
      filnam='SPCIN1'
      call molcas_binaryopen_vanilla(luTmp2,filnam)
c      Open(luTmp2,File=filnam,Form='unformatted',Status='unknown')
      ReWind (LuTmp1)
      ReWind (LuTmp2)
*
      Read (LuTmp1) nq,nQQ
      Write(LuTmp2) nq,nQQ
      Call GetMem('Temp_rK','Allo','Real',ipTmp,nQQ)
      Do iq = 1, nq
         Read (LuTmp1) qLbl,(Work(ipTmp+i),i=0,nQQ-1)
         Write(LuTmp2) qLbl,(Work(ipTmp+i),i=0,nQQ-1)
      End Do
      Call GetMem('Temp_rK','Free','Real',ipTmp,nQQ)
*
      Close  (LuTmp1)
      Close  (LuTmp2)
      End If
*
      Return
      End
