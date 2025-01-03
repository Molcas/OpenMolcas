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
      SubRoutine Ex_spin(rD,Fock,Temp1,ntemp,Temp2)
      use Constants, only: Zero,Half
      use MCLR_Data, only: nDens2, ipCM, nNA
      use input_mclr, only: nSym,nAsh,nIsh,nBas
      Implicit None
      Integer nTemp
      Real*8 rD(*),Fock(*),Temp1(nTemp),Temp2(nDens2)
      Integer jS, kS, llB, lB, jjB, jB
      Real*8 rDens
*
*
      Temp2(:)=Zero
      Do jS=1,nSym
         Call FZero(Fock(ipCm(js)),nbas(js)**2)
         Do kS=1,nsym
*
*...To be debugged!
*
            If (nBas(js)*nash(ks).gt.0) Then
*
               Do llB = 1, nAsh(kS)
                  lB=nIsh(kS)+llB
                  Do jjB = 1, nAsh(kS)
                     jB=nIsh(kS)+jjB
*
                     Call Exch(jS,kS,jS,kS,lB,jB,Temp2,Temp2)
                     rDens=-Half*rD(nna*(jjB-1)+llB)
                     Call DaXpY_(nBas(jS)**2,rDens,Temp1,1,
     &                          Fock(ipCM(jS)),1)

*
                  End Do
               End Do
*
            End If
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End SubRoutine Ex_spin
