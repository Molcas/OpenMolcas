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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine PLF_RICD(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                    iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                    TInt,nTInt,mTInt,iTOff,iOffij,iOffkl)
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
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "WrkSpc.fh"
*
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp), TInt(nTInt,mTInt)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4)
      Logical Shijij
      Common /ibas_ricd/ jbas_, lbas_
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
C     Call qEnter('PLF_RICD')
      irout = 109
      iprint = nprint(irout)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      r1=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
      r2=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
      Write (6,*) ' Sum=',r1
      Write (6,*) ' Dot=',r2
      Call RecPrt(' In PLF_RICD: AOInt',' ',
     &                              AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      iAOsti=iAOst(1)
      iAOstj=iAOst(2)
      iAOstk=iAOst(3)
      iAOstl=iAOst(4)
C     Write (6,*) 'iAOsti,iAOstj,iAOstk,iAOstl=',
C    &             iAOsti,iAOstj,iAOstk,iAOstl
      iAOi=iAO(1)
      iAOj=iAO(2)
      iAOk=iAO(3)
      iAOl=iAO(4)
C     Write (6,*) 'iAOs=',iAO
C     Write (6,*) 'kOps=',kOp
C     Write (6,*) 'iTOff,iOffij,iOffkl=',iTOff,iOffij,iOffkl
C     Write (*,*) 'iBas,jBas,kBas,lBas=',iBas,jBas,kBas,lBas
*
*     The writing of the integrals here are shell blocked.
*
      Do i1 = 1, iCmp
         iSOs(1)=iAOtSO(iAOi+i1,kOp(1))+iAOsti
         Do i2 = 1, jCmp
            iSOs(2)=iAOtSO(iAOj+i2,kOp(2))+iAOstj
            Do i3 = 1, kCmp
               iSOs(3)=iAOtSO(iAOk+i3,kOp(3))+iAOstk
               Do i4 = 1, lCmp
                  iSOs(4)=iAOtSO(iAOl+i4,kOp(4))+iAOstl
*
                iSO =iSOs(1)
                jSO =iSOs(2)
                kSO =iSOs(3)
                lSO =iSOs(4)
*
C               Write (6,*)
C               Write (6,*) 'i1,i2,i3,i4,iSOs=',i1,i2,i3,i4,iSOs
C               Write (6,*) 'iBas,jBas,kBas,lBas=',iBas,jBas,kBas,lBas
*
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   Do kSOk = kSO, kSO+kBas-1
                      If (iAO(3).eq.iAO(4)) Then
                         iSOkl=iTri(kSOk,lSOl) + iOffkl
                      Else
                         iSOkl=(kSOk-1)*lCmp*lBas_+ lSOl + iOffkl
                      End If
                      Do jSOj = jSO, jSO+jBas-1
                         Do iSOi = iSO, iSO+iBas-1
                            nijkl = nijkl + 1
                            AInt=AOint(nijkl,i1,i2,i3,i4)
                            If (iAO(1).eq.iAO(2)) Then
                               iSOij=iTri(iSOi,jSOj) + iOffij
                            Else
                               iSOij=(iSOi-1)*jCmp*jBas_+ jSOj + iOffij
                            End If
*
C                           Write (6,*) 'iSOij,iSOkl,AInt=',
C    &                                   iSOij,iSOkl,AInt
                            ijSOij=Max(iSOij,iSOkl)-iTOff
                            klSOkl=Min(iSOij,iSOkl)
                            TInt(klSOkl,ijSOij)= AInt
*
                         End Do
                      End Do
                   End Do
                End Do
*
               End Do
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('TInt','(45G8.2)',TInt,nTInt,mTInt)
#endif
*                                                                      *
************************************************************************
*                                                                      *
C     Call qExit('PLF_RICD')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
      End If
      End
