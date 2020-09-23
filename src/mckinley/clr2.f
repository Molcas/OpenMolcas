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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine Clr2(rIn,rOut,ibas,icmp,jbas,jcmp,
     &                iaoi,iaoj,naco,ishell,
     &                temp1,temp2,temp3,temp4,temp5,temp6)
*
      use pso_stuff
      use SOAO_Info, only: iAOtSO
      use Symmetry_Info, only: iOper
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "etwas.fh"
#include "info.fh"
#include "buffer.fh"
#include "disp.fh"
#include "disp2.fh"
#include "WrkSpc.fh"
*
      Real*8 rIn(ibas*icmp*jbas*jcmp,0:nIrrep-1,
     &       nAco*(1+naco)/2,*)
      real*8 rout(*)
      Real*8 Temp1(ibas,icmp,*)
      Real*8 Temp2(*)
      Real*8 Temp4(ibas,icmp,nACO)
      Real*8 Temp5(jbas,jcmp,nACO)
      Real*8 Temp3(jbas,jcmp,*),Temp6(*)
      integer ishell(4),na(0:7),ipp(0:7)
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

      call dcopy_(Naco**4,[0.0d0],0,Temp2,1)
      call dcopy_(nACO*ICMP*IBAS,[0.0d0],0,Temp4,1)
      call dcopy_(nACO*JCMP*JBAS,[0.0d0],0,Temp5,1)
      nnA=0
      Do iS=0,nIrrep-1
       nA(iS)=nNA
       nna=nna+nAsh(is)
      End Do
      n=0
      Do i=1,nirrep
       n=n+ldisp(i-1)
      end do
      n=ibas*icmp*jbas*jcmp*nIrrep*nAco*(1+naco)/2*n
*
      ni=iCmp*iBas
      nj=jCmp*jBas
      ipi=1


      ipj=ipi+naco*ibas*icmp

      Call PckMo2(temp6(ipi),nAcO,icmp,iBas,jcmp,jBas,iaoi,iaoj)
      id=0
      Do mIrr=0,nIrrep-1
       iiii=0
       Do iS=0,nIrrep-1
        js=nrOpr(ieor(ioper(is),iOper(mIrr)))
        ipp(is)=iiii
        iiii=nbas(is)*nash(js)+iiii
       End Do
       Do iDisp=1,lDisp(mIrr)
        iD=id+1
        ia=1
        Do iIrr=0,nIrrep-1
         kl=0
         k=0
         Do kIrr=0,nIrrep-1
          Do kAsh=1,nAsh(kIrr)
           k=k+1
           l=0
           Do lIrr=0,kIrr
            kls=iEOR(iOper(kIrr),iOper(lIrr))
            jIrr=nropr(ieor(iEOR(iOper(iIrr),iOper(mIrr)),kls))
            ja=1
            Do j=0,jirr-1
             ja=ja+nAsh(j)
            End Do
*
*           Symmetry of Q matrix
*
            iis=nropr(ieor(iOper(jIrr),ioper(mIrr)))
            jis=nropr(ieor(iOper(iIrr),ioper(mIrr)))
*
            lMax=nAsh(lIrr)
            If (lIrr.eq.kirr) lmax=kash
            Do lAsh=1,lMax
             l=lash+nA(lIrr)
             kl=itri(k,l)
*
*            id,iirr,jirr,kA,lA
*
             If (nash(jirr).ne.0)
     &       Call DGEMM_('N','N',
     &                   ni,nAsh(jIrr),nj,
     &                   1.0d0,rin(1,iIrr,kl,id),ni,
     &                   Temp6(ipj+(ja-1)*jcmp*jBas),nj,
     &                   0.0d0,Temp1,ni)
            If (nash(iirr).ne.0)
     &       Call DGEMM_('T','N',
     &                   nash(iIrr),nAsh(jIrr),ni,
     &                   1.0d0,Temp6(ipi+(ia-1)*icmp*ibas),ni,
     &                   Temp1,ni,
     &                   0.0d0,Temp2,nash(iirr))
*
*
             Do iC=1,iCmp
              Do iB=1,iBas
               Do i=1,nAsh(jis)
                ih=i+na(jis)
                Temp4(iB,ic,i)=0.0d0
                Do iAsh=1,nAsh(jirr)
                 jh=iash+na(jirr)
                 fact=1.0d00
                 iij=itri(ih,jh)
                 if(iij.ge.kl .and. k.eq.l) fact=2.0d00
                 if(iij.lt.kl .and. ih.eq.jh) fact=2.0d00
                  If (k.ne.l) FacT=fact*2.0d0
                 rd=G2(itri(iij,kl),1)*fact
                 Temp4(iB,ic,i)=Temp4(ib,ic,i)+
     &                 Temp1(ib,ic,iash)*rd
                End Do
               End Do
              End Do
             End Do
*
             ipF=ipDisp3(id)-1+ipp(iirr)
             Do jAsh=1,nAsh(jis)
              Do iC=1,iCmp
               lSO=iAOtSO(iAOi+iC,iIrr)
               If (lso.gt.0) Then
                Do iB=1,iBas
                 lsl=lSO+ib-1
                 ipFKL=ipF+(jAsh-1)*nBas(iIrr)+lsl
                 rout(ipFKL)=rout(ipFKL)+Temp4(ib,ic,jash)
                End Do
               End If
              End Do
             End DO
*
             If (iShell(1).ne.ishell(2)) Then
              If (nash(jirr).ne.0)
     &        Call DGEMM_('T','N',
     &                    nj,nAsh(jIrr),ni,
     &                    1.0d0,rin(1,jIrr,kl,id),ni,
     &                    Temp6(ipi+(ja-1)*icmp*ibas),ni,
     &                    0.0d0,Temp3,nj)
              If (nash(iirr).ne.0)

     &        Call DGEMM_('T','N',nAsh(iirr),nAsh(jirr),nj,
     &                  One,Temp6(ipj+(ia-1)*jcmp*jBas),nj,
     &                  Temp3,nj,
     &                  one,Temp2,nAsh(iirr))
*
              Do jC=1,jCmp
               Do jB=1,jBas
                Do i=1,nAsh(jis)
                 ih=i+na(jis)
                 Temp5(jB,jc,i)=0.0d0
                 Do iAsh=1,nAsh(jirr)
                  jh=iash+na(jirr)
                  fact=1.0d00
                  iij=itri(ih,jh)
                  if(iij.ge.kl .and. k.eq.l) fact=2.0d00
                  if(iij.lt.kl .and. ih.eq.jh) fact=2.0d00
                  If (k.ne.l) FacT=fact*2.0d0
                  rd=G2(itri(iij,kl),1)*fact
                  Temp5(jB,jc,i)=Temp5(jb,jc,i)+Temp3(jb,jc,iash)*rd
                 End Do
                End Do
               End Do
              End Do
*
              ipf=ipDisp3(id)-1+ipp(iirr)
              Do iAsh=1,nAsh(jis)
               Do jC=1,jCmp
                iSO=iAOtSO(iAOj+jC,iIrr)
                If (iso.gt.0) Then
                 Do jB=1,jBas
                  i=iSO+jb-1
                  ipFKL=ipF+(iAsh-1)*nBas(iIrr)+i
                  rout(ipFKL)=rout(ipFKL)+Temp5(jb,jc,iash)
                 End Do
                End If
               End Do
              End DO
             End If
*
* Distribute integrals
*
             if (iirr.ge.jirr) then
             klt=itri(k,l)
             Do iAsh=1,nash(iirr)
              ij1=iAsh+na(iirr)
              Do jAsh=1,nash(jirr)
               ij2=jAsh+na(jirr)
               If (ij1.ge.ij2) Then
               ij12=itri(ij1,ij2)
               If (ij12.le.klt) Then
               ipM=ipMO(iD,1)+iTri(ij12,klt)-1
               ipm2=iash+(jash-1)*nash(iirr)
               rOut(ipm)=rout(ipm)+Temp2(ipm2)
               End If
               End If
              End Do
             End Do
             end if
            End DO ! lash
           End DO ! lirr
          End DO ! kash
         End DO ! kirr
         ia=ia+nAsh(iIrr)
        End DO ! iirr
       End DO ! ndisp
      End DO ! msym
      Return
      End
