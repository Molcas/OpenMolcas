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
      SubRoutine Integral_WrOut_Cho(
     &                           iCmp,iShell,MapOrg,
     &                           iBas,jBas,kBas,lBas,kOp,
     &                           Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                           AOInt,SOInt,nSOint,
     &                           iSOSym,nSkal,nSOs,
     &                           TInt,nTInt,FacInt,itOffs,mSym,
     &                           Dens,Fock,LDens,ExFac,NDens,
     &                           ind,nind,FckNoClmb,FckNoExch)
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "cholesky.fh"

      Character*18 SecNam
      Parameter (SecNam = 'Integral_WrOut_Cho')
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        itOffs(0:mSym-1,0:mSym-1,0:mSym-1), MapOrg(4)
      Logical Shijij,IJeqKL,FckNoClmb,FckNoExch
      Real*8 Dens(lDens,nDens), Fock(lDens,nDens), ExFac(nDens)
      Integer Ind(nInd,nInd,2)
*
* some dummy assignments to avoid compiler warnings about unused
* variables.
*
      If (lDens.gt.0.and.nDens.gt.0.and.FckNoClmb.and.FckNoExch.and.
     &    nInd.gt.0.and.mSym.gt.0.and.nSkal.gt.0) Then
         xDummy_1  = Dens(1,1)
         xDummy_2  = Fock(1,1)
         xDummy_3  = FacInt
         xDummy_4  = ExFac(1)
         iDummy_1  = Ind(1,1,1)
         iDummy_2  = itOffs(0,0,0)
         iDummy_3  = MapOrg(1)
      End If
*
* call sorting routine
*
      If (IfcSew .eq. 1) Then
         If (Petite) Then
           Call PLF_Cho(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else If (IfcSew .eq. 2) Then
         If (Petite) Then
           Call PLF_Cho_2(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho_2(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else If (IfcSew .eq. 3) Then
         If (Petite) Then
           Call PLF_Cho_3(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho_3(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else
         Write(6,*)
         Write(6,*)
         Write(6,*) '!!!!!!!!!! IfcSew=',IfcSew,' !!!!!!!!!!'
         Call Cho_Quit('IfcSew out of bounds in '//SecNam,105)
      End If
*
      Return
      End
