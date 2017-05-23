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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
      subroutine handle2name(handle,name)
************************************************************
*
*   <DOC>
*     <Name>handle2name</Name>
*     <Syntax>Call handle2name(handle,name)</Syntax>
*     <Arguments>
*       \Argument{handle}{file handle}{int}{in}
*       \Argument{name}{file name}{char*}{out}
*     </Arguments>
*     <Purpose>retrive the file name from molcas I/O</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Valera Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The routine can be called from aixrd/wr
*        return the file name, or 'Unknown'
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Integer (a-z)
#include "switch.fh"
#include "ctl.fh"
      Character*(*) name
*----------------------------------------------------------------------*
       name='Unknown'
       do i=1,MxFile
        if(CtlBlk(pHndle,i).eq.handle) then
          name=FCtlBlk(i)
          return
        endif
       enddo
      End
