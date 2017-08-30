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
*  handle2name
*
*> @brief
*>   Retrieve the file name from molcas I/O
*> @author Valera Veryazov
*>
*> @details
*> The routine can be called from ::aixrd or ::aixwr.
*> Return the file name, or '``Unknown``'.
*>
*> @param[in]  handle file handle
*> @param[out] name   file name
************************************************************************
      subroutine handle2name(handle,name)
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
