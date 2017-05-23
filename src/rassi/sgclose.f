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
      Subroutine SGClose(iSGStruct)
#include "Struct.fh"
      Dimension iSGStruct(nSGSize)
C Unpack structure iSGStruct:
      nSym   =iSGStruct(1)
      nLev   =iSGStruct(2)
      lISm   =iSGStruct(3)
      nVert  =iSGStruct(4)
      lDRT   =iSGStruct(5)
      lDown  =iSGStruct(6)
      lUp    =iSGStruct(7)
      MidLev =iSGStruct(8)
      MVSta  =iSGStruct(9)
      MVEnd  =iSGStruct(10)
      lMAW   =iSGStruct(11)
      lLTV   =iSGStruct(12)
      Call GetMem('ISm','Free','Integer',lISm,nLev)
      Call GetMem('DRT','Free','Inte',lDRT,4*nVert)
      Call GetMem('Down','Free','Inte',lDown,4*nVert)
      Call GetMem('Up','Free','Inte',lUp,4*nVert)
      Call GetMem('MAW','Free','Inte',lMAW,4*nVert)
      Call GetMem('LTV','Free','Inte',lLTV,nLev+2)
      return
      end
