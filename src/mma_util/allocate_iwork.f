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
* Copyright (C) Roland Lindh                                           *
************************************************************************
      Subroutine Allocate_iWork(ip,n)
************************************************************
*
*   <DOC>
*     <Name>Allocate\_iWork</Name>
*     <Syntax>Call Allocate\_iWork(ip,n)</Syntax>
*     <Arguments>
*       \Argument{ip}{pointer to iWork}{integer}{out}
*       \Argument{n}{length}{integer}{in}
*     </Arguments>
*     <Purpose>To allocate integer memory of length n and
*              to return the corresponding pointer ip.</Purpose>
*     <Dependencies>GETMEM/Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     Short parameter list wrapper to GETMEM for allocation of
*     memory of type INTEGER. Label defaulted to "AiW".
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('AiW','Allo','Inte',ip,n)
c
c  Hey, if you want to debug a case, where allocate_iwork
c       takes memory, but never released it:
c  1. look at iPos printed by GetMem('LIST',,,,)
c  2. modify the numebr below
c  3. run with MOLCAS_DEBUGGER=gdb MOLCAS_BOMB=Yes
c
c      if(ip.eq.2878949) call abend
*
      Return
      End
