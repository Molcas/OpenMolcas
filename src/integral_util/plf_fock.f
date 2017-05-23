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
      Subroutine PLF_Fck1(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                     MapOrg,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,
     &                     kOp,D,F,ind,nind,ExFac,NoClmb,NoExch)
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
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp),
     &       D(*), F(*)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4), MapOrg(4),
     &        ind(nind,*)
      Logical Shijij, usShij, usShkl, NoExch, NoClmb
*
      usShij = iShell(1).eq.iShell(2)
      usShkl = iShell(3).eq.iShell(4)
      fac_c=One
      if(usShij) fac_c=fac_c*Half
      if(usShkl) fac_c=fac_c*Half
      if(Shijij) fac_c=fac_c*Half
      fac_e=fac_c*Exfac
      If (NoExch) fac_e=Zero
      If (NoClmb) fac_c=Zero
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
                     Do kSOk = kSO, kSO+kBas-1
                        kl=ind(kSOk,lSOl)
                        dkl=D(kl)*Four
                        fkl=Zero
                        Do jSOj = jSO, jSO+jBas-1
                           jk=ind(jSOj,kSOk)
                           jl=ind(jSOj,lSOl)
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
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(MapOrg)
      End
      Subroutine PLF_Fck2(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                     MapOrg,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,
     &                     kOp,D,F,ldens,ndens,ind,nind,ExFac,
     &                     NoClmb,NoExch)
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
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp),
     &       D(ldens,ndens), F(ldens,ndens),exfac(ndens)
      Integer iShell(4), iAO(4), kOp(4),
     &        iAOst(4), iSOs(4), MapOrg(4), ind(nind,*)
      Logical Shijij, usShij, usShkl, NoExch(ndens), NoClmb(ndens)
*
      usShij = iShell(1).eq.iShell(2)
      usShkl = iShell(3).eq.iShell(4)
      fac_c=One
      if(usShij) fac_c=fac_c*Half
      if(usShkl) fac_c=fac_c*Half
      if(Shijij) fac_c=fac_c*Half
      fac_e1=fac_c*Exfac(1)
      fac_e2=fac_c*Exfac(2)
      If (NoExch(1)) fac_e1=Zero
      If (NoExch(2)) fac_e2=Zero
      If (NoClmb(1)) fac_c=Zero
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
                     Do kSOk = kSO, kSO+kBas-1
                        kl=ind(kSOk,lSOl)
                        dkl=D(kl,1)*Four
                        fkl=Zero
                        Do jSOj = jSO, jSO+jBas-1
                           jk=ind(jSOj,kSOk)
                           jl=ind(jSOj,lSOl)
                           Do iSOi = iSO, iSO+iBas-1
                              nijkl = nijkl + 1
                              val=AOint(nijkl,i1,i2,i3,i4)
                              val_c=val*fac_c
                              val_e1=val*fac_e1
                              val_e2=val*fac_e2
                              ij=ind(iSOi,jSOj)
                              ik=ind(iSOi,kSOk)
                              il=ind(iSOi,lSOl)
                              fkl=fkl+val_c*D(ij,1)
                              F(ij,1)=F(ij,1)+val_c*dkl
                              F(ik,1)=F(ik,1)-val_e1*D(jl,1)
                              F(jl,1)=F(jl,1)-val_e1*D(ik,1)
                              F(il,1)=F(il,1)-val_e1*D(jk,1)
                              F(jk,1)=F(jk,1)-val_e1*D(il,1)
                              F(ik,2)=F(ik,2)-val_e2*D(jl,2)
                              F(jl,2)=F(jl,2)-val_e2*D(ik,2)
                              F(il,2)=F(il,2)-val_e2*D(jk,2)
                              F(jk,2)=F(jk,2)-val_e2*D(il,2)
                           End Do
                        End Do
                        F(kl,1)=F(kl,1)+fkl*Four
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(MapOrg)
      End
