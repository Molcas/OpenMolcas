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
      Subroutine Free_iWork(ip)
************************************************************
*
*   <DOC>
*     <Name>Free\_iWork</Name>
*     <Syntax>Call Free\_iWork(ip)</Syntax>
*     <Arguments>
*       \Argument{ip}{pointer to memory in iWork}i{integer}{inout}
*     </Arguments>
*     <Purpose>To deallocate memory in iWork associated with pointer ip.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Roland Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     The array associated to pointer/identified ip in iWork is deallocted.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('FiW','Free','Inte',ip,nDum)
*
      Return
      End
