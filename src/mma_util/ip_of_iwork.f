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
      Integer Function ip_of_iWork(A)
************************************************************
*
*   <DOC>
*     <Name>ip\_of\_iWork</Name>
*     <Syntax>ip\_of\_iWork(A)</Syntax>
*     <Arguments>
*       \Argument{A}{first element of integer array}{integer}{in}
*     </Arguments>
*     <Purpose>To return pointer to be used in iWork which correspond to
*              the first element of array A.</Purpose>
*     <Dependencies>iiLoc</Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     Function returns pointer such that A and iWork(ip\_of\_iWork(A)) refers to
*     the same memory location.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Integer A
*
      loc1=(iiLoc(A)-iiLoc(iWork(ip_iDummy)))
      loc2=(iiLoc(iWork(ip_iDummy+1))-iiLoc(iWork(ip_iDummy)))
      ip_of_iWork = ip_iDummy + loc1/loc2
*
      Return
      End
