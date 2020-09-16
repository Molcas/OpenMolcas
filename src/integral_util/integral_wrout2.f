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
      SubRoutine Integral_WrOut2(iCmp,iShell,MapOrg,
     &                           iBas,jBas,kBas,lBas,kOp,
     &                           Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                           AOInt,SOInt,nSOint,
     &                           iSOSym,nSkal,nSOs,
     &                           TInt,nTInt,FacInt,itOffs,nSym,
     &                           Dens,Fock,LDens,ExFac,NDens,
     &                           ind,nind,FckNoClmb,FckNoExch)
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        itOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4)
      Logical Shijij,IJeqKL,FckNoClmb,FckNoExch
*
      If (nSym.eq.1) Then
        Call PLF2(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &           iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &           iBas,jBas,kBas,lBas,kOp)
      Else
        Call IndSft2(iCmp,iShell,
     &               iBas,jBas,kBas,lBas,Shijij,
     &               iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(MapOrg)
         Call Unused_integer(nSkal)
         Call Unused_real_array(TInt)
         Call Unused_real(FacInt)
         Call Unused_integer_array(itOffs)
         Call Unused_real(Dens)
         Call Unused_real(ExFac)
         Call Unused_real(Fock)
         Call Unused_integer(LDens)
         Call Unused_integer(NDens)
         Call Unused_integer(ind)
         Call Unused_integer(nind)
         Call Unused_logical(FckNoClmb)
         Call Unused_logical(FckNoExch)
      End If
      End
