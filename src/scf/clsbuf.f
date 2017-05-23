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
* Copyright (C) 1998, Roland Lindh                                     *
************************************************************************
      SubRoutine ClsBuf(nDisc,nCore)
************************************************************************
*                                                                      *
*  Object: Close I/O buffer fro semi-direct SCF                        *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, Sweden. October '98                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "IOBuf.fh"
*
      If (OnDisk) Call EAFClose(LuTmp)
      If (nCore.ne.0) Call GetMem('IOBuf','Free','Real',ipBuf,lBuf*nBuf)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nDisc)
      End
