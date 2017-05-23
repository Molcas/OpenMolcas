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
      Subroutine Get_AnalHess(ipAnalHess,nAnalHess)
************************************************************
*
*   <DOC>
*     <Name>Get\_AnalHess</Name>
*     <Syntax>Call Get\_AnalHess(ipAnalHess,nAnalHess)</Syntax>
*     <Arguments>
*       \Argument{ipAnalHess}{pointer to array with the symmetry blocked nuclear Hessian
*                             in cartesian coordinates}{Integer}{out}
*       \Argument{nAnalHess}{size of the array of the symmetry blocked nuclear Hessian}{Integer}{out}
*     </Arguments>
*     <Purpose>To read the the symmetry blocked nuclear Hessian from the run file and to return a
*              pointer to the arrays location in Work and the length of the array.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The utility will the the symmetry blocked nuclear Hessian from the run file and
*                  return a pointer and the length of the array.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found

      Label='Analytic Hessian'
      Call qpg_dArray(Label,Found,nAnalHess)
      If(.not.Found .or. nAnalHess.eq.0) Return
      Call GetMem('AnalHess','Allo','Real',ipAnalHess,nAnalHess)
      Call Get_dArray(Label,Work(ipAnalHess),nAnalHess)

      Return
      End
