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
      Subroutine Allocate_Work(ip,n)
************************************************************
*
*   <DOC>
*     <Name>Allocate\_Work</Name>
*     <Syntax>Call Allocate\_Work(ip,n)</Syntax>
*     <Arguments>
*       \Argument{ip}{pointer to Work}{integer}{out}
*       \Argument{n}{length}{integer}{in}
*     </Arguments>
*     <Purpose>To allocate real memory of length n and
*              to return the corresponding pointer ip.</Purpose>
*     <Dependencies>GETMEM</Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     Short parameter list wrapper to GETMEM for allocation of
*     memory of type REAL. Label defaulted to "AW".
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('AW','Allo','Real',ip,n)
*
      Return
      End
