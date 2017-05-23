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
      Subroutine FetchTDM(nB,nS,iBigT,TDMchar)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension iTocBig(MxStOT)

      Character TDMchar*6

      iDisk=0
      kaunter=0
      nSize=nB*(nB+1)/2
      index=0
      Lu=72
      Lu=IsFreeUnit(Lu)
      Call DaName(Lu,TDMchar)
      Call iDaFile(Lu,2,iTocBig,MxStOT,iDisk)
      Do 99991, iS1=1,nS
        Do 99992, iS2=1,iS1
          kaunter=kaunter+1
          iDisk=iTocBig(kaunter)
          Call dDaFile(Lu,2,Work(iBigT+index),nSize,iDisk)
          index=index+nSize
99992   Continue
99991 Continue
      Call DaClos(Lu)

      Return
      End
