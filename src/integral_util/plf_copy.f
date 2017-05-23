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
* Copyright (C) 1990,1998, Roland Lindh                                *
************************************************************************
      Subroutine PLF_Copy(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                    MapOrg,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,
     &                    kOp,TInt,nTint,FacInt,
     &                    nShi,nShj,nShk,nShl,
     &                    nShOffi,nShOffj,nShOffk,nShOffl,
     &                    PL_BAddr_Inc)
************************************************************************
*                                                                      *
*  object: to shift and index the petite list format integrals.        *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          May '90                                                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
*
      External PL_BAddr_Inc
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp), TInt(nTInt)
      Integer iShell(4), iAO(4), kOp(4),
     &        iAOst(4), iSOs(4), MapOrg(4)
      Logical Shijij
*
*     Call qEnter('PLF_Copy')
      irout = 109
      iprint = nprint(irout)
      If (iPrint.ge.49) Then
         r1=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
         r2=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
         Write (6,*) ' Sum=',r1
         Write (6,*) ' Dot=',r2
      End If
      If (iPrint.ge.99) Call RecPrt(' In Plf_Copy: AOInt',' ',
     &                              AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
*
      iAOsti=iAOst(1)
      iAOstj=iAOst(2)
      iAOstk=iAOst(3)
      iAOstl=iAOst(4)
      iAOi=iAO(1)
      iAOj=iAO(2)
      iAOk=iAO(3)
      iAOl=iAO(4)
*
      Call PL_BAddr_Inc(MapOrg,nShi,nShj,nShk,nShl,
     &                  nShOffi,nShOffj,nShOffk,nShOffl,
     &                  Inci,Incj,Inck,Incl,iOff4)
*
      Do i4 = 1, lCmp
         iSOs(4)=iAOtSO(iAOl+i4,kOp(4))+iAOstl
         Do i3 = 1, kCmp
            iSOs(3)=iAOtSO(iAOk+i3,kOp(3))+iAOstk
            Do i2 = 1, jCmp
               iSOs(2)=iAOtSO(iAOj+i2,kOp(2))+iAOstj
               Do i1 = 1, iCmp
                  iSOs(1)=iAOtSO(iAOi+i1,kOp(1))+iAOsti
*
                  iSO =iSOs(1)
                  jSO =iSOs(2)
                  kSO =iSOs(3)
                  lSO =iSOs(4)
*
                  nijkl = 0
                  Do lSOl = lSO, lSO+lBas-1
                     iAdr4l = iOff4 + lSOl*Incl
                     Do kSOk = kSO, kSO+kBas-1
                        iAdr4kl = iAdr4l + kSOk*Inck
                        Do jSOj = jSO, jSO+jBas-1
                           iAdr4jkl = iAdr4kl + jSOj*Incj
                           Do iSOi = iSO, iSO+iBas-1
                              iAdr4ijkl = iAdr4jkl + iSOi*Inci
                              nijkl = nijkl + 1
                              TInt(iAdr4ijkl)
     >                        =AOint(nijkl,i1,i2,i3,i4)*FacInt
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
*
*     Call qExit('PLF_Copy')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
      End If
      End


      SubRoutine PL_BAddr_Inc_ijkl(MapOrg,
     &                             nShi,nShj,nShk,nShl,
     &                             nShOffi,nShOffj,nShOffk,nShOffl,
     &                             Inci,Incj,Inck,Incl,iBAddr)
*     compute addresses for integrals in order ijkl...
      Implicit Real*8 (A-H,O-Z)
*     declaration of subr parameters...
      Integer MapOrg(4),
     &        nShi,nShj,nShk,nShl,nShOffi,nShOffj,nShOffk,nShOffl,
     &        Inci,Incj,Inck,Incl,iBAddr
*     declaration of local variables...
      Integer Inc(4),IncP(4)
*
      Inc(4)=1
      Inc(3)=nShl
      Inc(2)=nShl*nShk
      Inc(1)=nShl*nShk*nShj
      iBAddr = 1 -Inc(4)*(nShOffl+1)-Inc(3)*(nShOffk+1)
     &           -Inc(2)*(nShOffj+1)-Inc(1)*(nShOffi+1)
      Do i = 1, 4
         IncP(i)=Inc(MapOrg(i))
      End Do
      Inci = IncP(1)
      Incj = IncP(2)
      Inck = IncP(3)
      Incl = IncP(4)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nShi)
      End


      SubRoutine PL_BAddr_Inc_jikl(MapOrg,
     &                             nShi,nShj,nShk,nShl,
     &                             nShOffi,nShOffj,nShOffk,nShOffl,
     &                             Inci,Incj,Inck,Incl,iBAddr)
*     compute addresses for integrals in order jikl...
      Implicit Real*8 (A-H,O-Z)
*     declaration of subr parameters...
      Integer MapOrg(4),
     &        nShi,nShj,nShk,nShl,nShOffi,nShOffj,nShOffk,nShOffl,
     &        Inci,Incj,Inck,Incl,iBAddr
*     declaration of local variables...
      Integer Inc(4),IncP(4)
*
      Inc(4)=1
      Inc(3)=nShl
      Inc(2)=nShl*nShk*nShi
      Inc(1)=nShl*nShk
      iBAddr = 1 -Inc(4)*(nShOffl+1)-Inc(3)*(nShOffk+1)
     &           -Inc(2)*(nShOffj+1)-Inc(1)*(nShOffi+1)
      Do i = 1, 4
         IncP(i)=Inc(MapOrg(i))
      End Do
      Inci = IncP(1)
      Incj = IncP(2)
      Inck = IncP(3)
      Incl = IncP(4)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nShj)
      End
