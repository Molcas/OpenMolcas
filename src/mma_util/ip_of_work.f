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
      Integer Function ip_of_work(A)
************************************************************
*
*   <DOC>
*     <Name>ip\_of\_Work</Name>
*     <Syntax>ip\_of\_Work(A)</Syntax>
*     <Arguments>
*       \Argument{A}{first element of real array}{real}{in}
*     </Arguments>
*     <Purpose>To return pointer to be used in Work which correspond to
*              the first element of array A.</Purpose>
*     <Dependencies>idLoc</Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     Function returns pointer such that A and Work(ip\_of\_Work(A)) refers to
*     the same memory location.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 A
*
      loc1=(idLoc(A)-idLoc(Work(ip_Dummy)))
      loc2=(idLoc(Work(ip_Dummy+1))-idLoc(Work(ip_Dummy)))
      ip_of_Work = ip_Dummy + loc1/loc2
*
      Return
      End
