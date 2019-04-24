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
      SubRoutine MakeMO(AOInt,Temp,nTemp,nInt,
     &                  MOInt,nMOInt,
     &                  iCmp,iCmpa,
     &                  ibasi,jbasj,kbask,lbasl,
     &                  nGr,index,
     &                  moip,naco,nop,indgrd,
     &                  ishll,ishell,rmoin,nmoin,iuvwx,iao,iaost,
     &                  buffer,ianga,c)
*
*   this is the driver for the two index transformation
*   it is not very efficent, but on the other hand it
*   usually to take more than a few percent of the total
*   CPU time, if someone notice something else I will
*   rewrite it, in the mean time, dont worry.
*


      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "buffer.fh"
*
*
      Logical pert(0:7),lc
      Integer iCmpa(4),
     &         index(3,4),ipPert(0:7),icmp(4),ibas(4),
     &         indgrd2(3,4,0:7),indgrd(3,4,0:nirrep-1),
     &         moip(0:nIrrep-1),nop(4),ishell(4),iuvwx(4),
     &         iao(4),iAOST(4),ianga(4),ishll(4)
      Real*8 Temp(nTemp),AOInt(nInt),rmoin(nmoin),MOInt(nMOInt),
     &       C(12),buffer(*)
*
      iMax=0
      nMax=0
      mSum=0
      nabcd=iBasi*jBasj*kBask*lBasl
      nijkl=icmp(1)*icmp(2)*icmp(3)*icmp(4)
      iBas(1)=iBasi
      iBas(2)=jBasj
      iBas(3)=kBask
      iBas(4)=lBasl
      Do ii=1,4
       imax=Max(iBas(ii)*iCmp(ii),imax)
       mSum=mSum+iBas(ii)*iCmp(ii)
      End Do
      imax=Max(iMax,nAco)
*
      nInt2=nabcd*nijkl

      ip=1
      ip0=ip
      ip=ip+nGr*nijkl*nabcd
      ip1=ip
      nScrtch=imax**4
      ip=ip+nScrtch
      ip2=ip
      ip=ip+nScrtch
      ip3=ip
      ip=ip+nScrtch
      ip5=ip
      ip=ip+iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)*
     &   iBas(1)*iBas(2)*iBas(3)*iBas(4)
      If (ip-1.gt.nTemp) Then
         Write (6,*) 'MakeMO: ip-1.gt.nTemp'
         Write (6,*) 'ip,nTemp=',ip,nTemp
         Call QTrace
         Call Abend()
      End If
*     ip=2
*     Temp(ip-1)=0.0d0
*     ip0=ip
*     ip=ip+nGr*nijkl*nabcd+1
*     Temp(ip-1)=0.0d0
*     ip1=ip
*     nScrtch=imax**4+1
*     ip=ip+nScrtch
*     Temp(ip-1)=0.0d0
*     ip2=ip
*     ip=ip+nScrtch
*     Temp(ip-1)=0.0d0
*     ip3=ip
*     ip=ip+nScrtch
*     Temp(ip-1)=0.0d0
*     ip5=ip
*     ip=ip+iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)*
*    &   iBas(1)*iBas(2)*iBas(3)*iBas(4)+1
*     Temp(ip-1)=0.0d0
*
      ipc=1
      ipD=ipc
*     nD=nACO**4
*     ipC=ipC+nd
      ipci=ipc
*     nCi=iBas(1)*iCmp(1)*nACO
      ipcj=ipc
*     nCj=iBas(2)*iCmp(2)*nACO
      ipck=ipc
      nCk=iBas(3)*iCmp(3)*nACO
      ipc=ipc+nCk
      ipcl=ipc
      nCl=iBas(4)*iCmp(4)*nAcO
      ipc=ipc+nCl
      nij=iCmp(1)*iBas(1)*iBas(2)*iCmp(2)
      nkl=iCmp(3)*iBas(3)*iBas(4)*iCmp(4)
      If (ipc-1.ne.nMoIn) Then
         Write (6,*) 'MakeMO: ipc-1.ne.nMoIn'
         Write (6,*) 'ipc,nMoIn=',ipc,nMoIn
         Call QTrace
         Call Abend()
      End If

*
      Call Sort_mck(AOInt,Temp(ip0),
     &          iBas(1),iBas(2),iBas(3),iBas(4),
     &          iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &          iBas(1),iBas(2),iBas(3),iBas(4),
     &          iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),
     &          nGr,nop,ianga,
     &          indgrd,indgrd2,ishll,C)
*
      Do iCent=1,4
       lc=.false.
       Do iCar=1,3
        Call ICopy(nIrrep,[0],0,ipPert,1)
        Call lCopy(nIrrep,[.false.],0,pert,1)
        lc=.false.
        Do iIrrep=0,nIrrep-1
         If (indgrd(icar,icent,iIrrep).ne.0) Then
            ipPert(iIrrep)=ipMO(abs(indgrd(iCar,iCent,iIrrep)),1)
            pert(iIrrep)=.true.
         End If
        If (IndGrd(iCar,icent,iirrep).ne.0) lC=.true.
        End Do
        If (lc) Then
        If (Index(iCar,iCent).gt.0) Then

*
         iGr=index(icar,icent)
         Call MOAcc(Temp(ip0+(iGr-1)*nijkl*nabcd),nINT2,
     &           Temp(ip1),Temp(ip2),Temp(ip3),nScrtch,
     &           MOInt,nMOINt,ishell,
     &           rmoin(ipCi),nCi,rmoin(ipCj),nCj,
     &           rmoin(ipCk),nCk,rmoin(ipCl),nCl,
     &           Moip,nACO,pert,nOp,ibas,icmpa,
     &           iCar,icent,indgrd,rmoin(ipD),
     &           DBLE(iuvwx(iCent))/DBLE(nIrrep),iao,iaost,
     &           buffer,Temp(ip2),nij,nkl,
     &           nBasis(ishll(1)),nBasis(ishll(2)),icmpa(1),icmpa(2))
*
        Else If (Index(iCar,iCent).lt.0) Then
         call dcopy_(nabcd*nijkl,[Zero],0,Temp(ip5),1)
         Do iCnt=1,4
           iGr=Index(iCar,iCnt)
           If (iGr.gt.0)Then
             ipFin=(iGr-1)*nijkl*nabcd+ip0
             Do ii=1,nabcd*nijkl
                Temp(ip5+ii-1)=Temp(ip5+ii-1)-Temp(ipFin-1+ii)
             End Do
           End If
         End Do
         Call MOAcc(Temp(ip5),nInt2,
     &           Temp(ip1),Temp(ip2),Temp(ip3),nScrtch,
     &           MOInt,nMOINt,ishell,
     &           rmoin(ipCi),nCi,rmoin(ipCj),nCj,
     &           rmoin(ipCk),nCk,rmoin(ipCl),nCl,
     &           moip,nACO,pert,nOp,ibas,icmpa,
     &           iCar,icent,indgrd,rMoin(ipD),
     &           DBLE(iuvwx(iCent))/DBLE(nIrrep),iao,iaost,
     &           buffer,Temp(ip2),nij,nkl,
     &           nBasis(ishll(1)),nBasis(ishll(2)),icmpa(1),icmpa(2))
*
        End If
        End If
       End Do
      End Do
      Return
      End
