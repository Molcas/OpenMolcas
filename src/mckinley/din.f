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
      use pso_stuff
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
#include "WrkSpc.fh"
      Real*8 Dens(nDens)
*                                                                      *
************************************************************************
*                                                                      *
      nTemp2=0
      Do iIrr=0,nIrrep-1
         nTemp2=Max(nTemp2,nBas(iIrr))
      End Do
      Call GetMem('Temp2','ALLO','REAL',ipTemp2,nTemp2**2)
      ip=ipCMO
      ipD=0
      Do iIrr=0,nIrrep-1
         If (nBas(iIrr).ne.0) Then
            Call DGEMM_('N','T',
     &                  nBas(iIrr),nBas(iIrr),nIsh(iIrr),
     &                  1.0d0,Work(ip),nBas(iIrr),
     &                  Work(ip),nBas(iIrr),
     &                  0.0d0,Work(ipTemp2),nBas(iIrr))
            Do iBas=1,nBas(iIrr)
               Do jBas=1,iBas-1
                  ip1=(iBas-1)*nBas(iIrr)+jBas
                  ip2=iBas*(iBas-1)/2+jBas
                  Dens(ipD+ip2)=Work(ipTemp2+ip1-1)*Four
               End Do
               ip1=(iBas-1)*nBas(iIrr)+iBas
               ip2=iBas*(iBas+1)/2
               Dens(ipD+ip2)=Work(ipTemp2+ip1-1)*Two
            End Do
            ip=ip+nBas(iIrr)**2
            ipd=ipD+nBas(iIrr)*(nBas(iIrr)+1)/2
         End If
      End Do
      Call GetMem('Temp2','FREE','REAL',ipTemp2,nTemp2**2)
      Return
      End
