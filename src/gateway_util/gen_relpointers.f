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
      Subroutine Gen_RelPointers(ibase)
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
*
      Do i = 1, MxShll
         ipCff(i)        = ipCff(i)         + ibase
         ipCff_Cntrct(i) = ipCff_Cntrct(i)  + ibase
         ipCff_Prim(i)   = ipCff_Prim(i)    + ibase
         ipExp(i)        = ipExp(i)         + ibase
         ipBk(i)         = ipBk(i)          + ibase
         ip_Occ(i)       = ip_Occ(i)        + ibase
         ipAkl(i)        = ipAkl(i)         + ibase
         ipFockOp(i)     = ipFockOp(i)      + ibase
      End Do
      Do i = 1, Mxdbsc
         ipCntr(i)     = ipCntr(i)     + ibase
         ipM1xp(i)     = ipM1xp(i)     + ibase
         ipM2xp(i)     = ipM2xp(i)     + ibase
         ipM1cf(i)     = ipM1cf(i)     + ibase
         ipM2cf(i)     = ipM2cf(i)     + ibase
         ipFragEner(i) = ipFragEner(i) + ibase
         ipFragCoef(i) = ipFragCoef(i) + ibase
         ipFragCoor(i) = ipFragCoor(i) + ibase
         ipPAM2xp(i)   = ipPAM2xp(i)   + ibase
         ipPAM2cf(i)   = ipPAM2cf(i)   + ibase
      End Do
      ipEF = ipEF + ibase
      ipOAM= ipOAM+ ibase
      ipOMQ= ipOMQ+ ibase
      ipDMS= ipDMS+ ibase
      ipWel= ipWel+ ibase
      ipAMP= ipAMP+ ibase
      ipRP1=ipRP1 + ibase
      ipXF = ipXF + ibase
      ipXMolnr=ipXMolnr + RtoI*ibase
      ipXEle=ipXEle + RtoI*ibase
*
      If (Allocated(iSD)) Then
         Call Nr_Shells(nSkal)
         Do iSkal = 1, nSkal
            iSD(4,iSkal)= iSD(4,iSkal) + ibase
            iSD(6,iSkal)= iSD(6,iSkal) + ibase
            iSD(8,iSkal)= iSD(8,iSkal) + ibase
         End Do
      End If
*
      Return
      End
