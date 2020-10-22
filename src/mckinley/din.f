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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine DIN(Dens)
*
      use Basis_Info, only: nBas
      use pso_stuff
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "etwas.fh"
#include "stdalloc.fh"
      Real*8 Dens(nDens)
      Real*8, Allocatable:: Temp2(:)
*                                                                      *
************************************************************************
*                                                                      *
      nTemp2=0
      Do iIrr=0,nIrrep-1
         nTemp2=Max(nTemp2,nBas(iIrr))
      End Do

      Call mma_allocate(Temp2,nTemp2**2,Label='Temp2')

      ip=1
      ipD=0
      Do iIrr=0,nIrrep-1

         If (nBas(iIrr)==0) Cycle

         Call DGEMM_('N','T',
     &               nBas(iIrr),nBas(iIrr),nIsh(iIrr),
     &               One,CMO(ip,1),nBas(iIrr),
     &                     CMO(ip,1),nBas(iIrr),
     &               Zero,Temp2,nBas(iIrr))
         Do iBas=1,nBas(iIrr)
            Do jBas=1,iBas-1
               ip1=(iBas-1)*nBas(iIrr)+jBas
               ip2=iBas*(iBas-1)/2+jBas
               Dens(ipD+ip2)=Temp2(ip1)*Four
            End Do
            ip1=(iBas-1)*nBas(iIrr)+iBas
            ip2=iBas*(iBas+1)/2
            Dens(ipD+ip2)=Temp2(ip1)*Two
         End Do
         ip=ip+nBas(iIrr)**2
         ipd=ipD+nBas(iIrr)*(nBas(iIrr)+1)/2

      End Do

      Call mma_deallocate(Temp2)

      Return
      End
