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
      Subroutine PickPoints(nPick,ipPick,ipDPick,nEPP,ipEPCo,Coo
     &                     ,dLimmo,BS)
      Implicit real*8(a-h,o-z)

#include "WrkSpc.fh"

      Dimension Coo(3),dLimmo(2)

      nPick=0
      Do iP=1,nEPP
        xtwo=(Work(ipEPCo+(iP-1)*3+0)-Coo(1))**2
        ytwo=(Work(ipEPCo+(iP-1)*3+1)-Coo(2))**2
        ztwo=(Work(ipEPCo+(iP-1)*3+2)-Coo(3))**2
        Distad=sqrt(xtwo+ytwo+ztwo)
        If(Distad.lt.dLimmo(2)*BS.and.Distad.gt.dLimmo(1)*BS) then
          iWork(ipPick+nPick)=iP
          Work(ipDPick+nPick)=Distad
          nPick=nPick+1
        Endif
      Enddo

      Return
      End
