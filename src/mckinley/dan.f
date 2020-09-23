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
      SubRoutine DAN(Dens)
*
      use Basis_Info, only: nBas
      use pso_stuff
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
#include "WrkSpc.fh"
      Real*8 Dens(nDens)
      Integer na(0:7),ipcm(0:7)
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      ipD=0
      nnA=0
      ndenssq=0
      ipCC=1
      Do i=0,nIrrep-1
         nDenssq=ndenssq+nBas(i)**2
         nA(i)=nnA
         ipcm(i)=ipCC
         nnA=nnA+nAsh(i)
         ipCC=ipCC+nBAs(i)**2
      End Do
*
      Call GetMem('Temp1','ALLO','REAL',ipTemp1,ndensSQ)
      Call GetMem('Temp2','ALLO','REAL',ipTemp2,ndensSQ)
      Call GetMem('Temp3','ALLO','REAL',ipUrk,ndensSQ)
*
      Do iS=0,nIrrep-1
         call dcopy_(nBas(is)**2,[0.0d0],0,Work(ipTemp1),1)
         If (nBas(is).gt.0) Then
            Do iB=1,nAsh(iS)
               iiB=nA(iS)+iB
               Do jB=1,nAsh(iS)
                  jjB=nA(iS)+jB
                  ijB=iTri(iiB,jjB)
                  ip1=ipTemp1+
     &                nBas(iS)*(nISh(iS)+iB-1)+nIsh(is)+jb -1
                  Work(ip1)=G1(ijB,1)
               End Do
            End Do
*
            Call DGEMM_('N','N',
     &                  nBas(is),nBas(is),nBas(is),
     &                  1.0d0,CMO(ipCM(iS),1),nBas(is),
     &                  Work(ipTemp1),nBas(is),
     &                  0.0d0,Work(ipUrk),nBas(is))
            Call DGEMM_('N','T',
     &                  nBas(is),nBas(is),nBas(is),
     &                  1.0d0,Work(ipUrk),nBas(is),
     &                  CMO(ipCM(is),1),nBas(is),
     &                  0.0d0,Work(ipTemp2),nBas(is))
*
            Do iBas=1,nBas(iS)
               Do jBas=1,iBas
                  ip1=(iBas-1)*nBas(iS)+jBas
                  ip2=iTri(iBas,jBas)
                  Fact=2.0d0
                  If (iBas.eq.jbas) Fact=1.0d0
                  Dens(ipD+ip2)=Work(ipTemp2+ip1-1)*Fact
               End Do
            End Do
            ipD=ipD+nBas(is)*(nBAs(is)+1)/2
         End If
      End Do
      Call GetMem('Temp1','FREE','REAL',ipTemp1,ndensSQ)
      Call GetMem('Temp2','FREE','REAL',ipTemp2,ndensSQ)
      Call GetMem('Temp3','FREE','REAL',ipUrk,ndensSQ)
*
      Return
      End
