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
      Subroutine PLF_CpFck(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                     MapOrg,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,
     &                     kOp,TInt,nTint,FacInt,D,F,ind,nind,NoExch,
     &                     nShi,nShj,nShk,nShl,
     &                     nShOffi,nShOffj,nShOffk,nShOffl,
     &                     PL_BAddr_Inc)
************************************************************************
*                                                                      *
*  object: to sift and index the petite list format integrals.         *
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
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp), TInt(nTInt),
     &       D(*), F(*)
      Integer iShell(4), iAO(4), kOp(4),
     &        iAOst(4), iSOs(4), MapOrg(4), ind(nind,*)
      Logical Shijij, usShij, usShkl, NoExch
*
*     Call qEnter('PLF_CpFck')
      irout = 109
      iprint = nprint(irout)
c     If (iPrint.ge.49) Then
c        r1=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
c        r2=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
c        Write (*,*) ' Sum=',r1
c        Write (*,*) ' Dot=',r2
c     End If
c     If (iPrint.ge.99) Call RecPrt(' In Plf_Copy: AOInt',' ',
c    &                              AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
*
*     and set up corresponding logicals...
      usShij = iShell(1).eq.iShell(2)
      usShkl = iShell(3).eq.iShell(4)
      fac_c=One
      if(usShij) fac_c=fac_c*Half
      if(usShkl) fac_c=fac_c*Half
      if(Shijij) fac_c=fac_c*Half
      fac_e=fac_c
      If (NoExch) fac_e=Zero
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
                        kl=ind(kSOk,lSOl)
                        dkl=D(kl)*Four
                        fkl=Zero
                        iAdr4kl = iAdr4l + kSOk*Inck
                        Do jSOj = jSO, jSO+jBas-1
                           jk=ind(jSOj,kSOk)
                           jl=ind(jSOj,lSOl)
                           iAdr4jkl = iAdr4kl + jSOj*Incj
                           Do iSOi = iSO, iSO+iBas-1
                              nijkl = nijkl + 1
                              val=AOint(nijkl,i1,i2,i3,i4)
                              val_c=val*fac_c
                              val_e=val*fac_e
                              ij=ind(iSOi,jSOj)
                              fkl=fkl+val_c*D(ij)
                              F(ij)=F(ij)+val_c*dkl
                              ik=ind(iSOi,kSOk)
                              il=ind(iSOi,lSOl)
                              F(ik)=F(ik)-val_e*D(jl)
                              F(jl)=F(jl)-val_e*D(ik)
                              F(il)=F(il)-val_e*D(jk)
                              F(jk)=F(jk)-val_e*D(il)
                              iAdr4ijkl = iAdr4jkl + iSOi*Inci
                              TInt(iAdr4ijkl)=val*FacInt
                           End Do
                        End Do
                        F(kl)=F(kl)+fkl*Four
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
*
*     Call qExit('PLF_CpFck')
      Return
      End
