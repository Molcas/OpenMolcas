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
      use Basis_Info
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
         ipPAM2xp(i)   = ipPAM2xp(i)   + ibase
         ipPAM2cf(i)   = ipPAM2cf(i)   + ibase
      End Do
*
      If (Allocated(iSD)) Then
         Call Nr_Shells(nSkal)
         Do iSkal = 1, nSkal
            iSD(4,iSkal)= iSD(4,iSkal) + ibase
            iSD(6,iSkal)= iSD(6,iSkal) + ibase
         End Do
      End If
*
      Return
      End
