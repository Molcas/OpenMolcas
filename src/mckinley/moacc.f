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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
      SubRoutine MOAcc(AOInt,nint,Temp1,Temp2,Temp3,nTemp,
     &      MOInt,nMO,ishell,
     &      Ci,nCi,Cj,nCj,Ck,nCk,Cl,nCl,moip,nACO,pert,nOp,
     &      iBasa,iCmpa,icar,icnt,indgrd,D,fact,iao,iaost,
     &      Buffer,Tempi,nij,nkl,nbasi,nbasj,icmp,jcmp)
************************************************************************
*                                                                      *
*     Transforms a batch of unsymmetrized integrals to                 *
*     active integral batches and FM                                   *
*     All MO combinations are constructed                              *
*     They will be needed with unsymmetric perurbations                *
*                                                                      *
*     Author: Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
*             University of Lund, Sweden. Januar '96                   *
************************************************************************
      use Symmetry_Info, only: iChTbl
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "etwas.fh"
#include "info.fh"
c#include "print.fh"
      Real*8 AOInt(nkl,nij),MOint(nMO),
     &       Temp1(nTemp),Temp2(naco,naco),
     &       Ck(nCk),Cl(nCl),D(*),
     &       Buffer(nbasi,icmp,nbasj,jcmp,0:nirrep-1,
     &       nAco*(naco+1)/2,*)
      Integer moip(0:nIrrep-1),nOp(4),
     &          ishell(4),iao(4),iAOST(4),
     &          ibasa(4),icmpa(4),indgrd(3,4,0:nirrep-1)
      Logical pert(0:nIrrep-1)
      Real*8 Prmt(0:7)

      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/

*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iCB=2**(icar-1)
      rFact=xPrmt(ioper(nOp(icnt)),icb)*fact
*
      iBas=iBasa(1)
      jBas=iBasa(2)
      kBas=iBasa(3)
      lBas=iBasa(4)
*
      kCmp=iCmpa(3)
      lCmp=iCmpa(4)
      kk=0
      Do kIrrep=0,nIrrep-1
       sfact=DBLE(ichtbl(kirrep,nop(3)))
       Do kAsh=1,nAsh(kIrrep)
        Do k=1,kcmp*kbas
         kk=kk+1
         Ck(kk)=Ck(kk)*sFact
        End Do
       End Do
      End Do
      kk=0
      Do kIrrep=0,nIrrep-1
       sfact=DBLE(ichtbl(kirrep,nop(4)))
       Do kAsh=1,nAsh(kIrrep)
        Do k=1,lcmp*lbas
         kk=kk+1
         Cl(kk)=Cl(kk)*sFact
        End Do
       End Do
      End Do
      rk = DNrm2_(nck,ck,1)
      rl = DNrm2_(ncl,cl,1)
      ij=0
      nt=lBas*lCmp*kcmp*kbas
      Do jc=1,jcmp
       Do jb=1+iaost(2),iaost(2)+jbas
        Do ic=1,icmp
        Do ib=1+iaost(1),iaost(1)+ibas
*
          ij=ij+1
          vij = DNrm2_(nt,AOInt(1,ij),1)
          If (Abs(vij*rk*rl).lt.cutint) goto 1000
          ipC=0
          Do kAsh=1,nAco
           ipM=(kAsh-1)*lbas*lcmp
           il=0
           Do i=1,lbas*lcmp
            Temp1(i+ipM)=0.0d0
            Do k=1,kCmp*kBas
             il=il+1
             Temp1(i+ipM)=Temp1(ipm+i)+Ck(k+ipc)*AOINT(il,ij)
            End Do
           End Do
           ipC=ipC+kBas*kCmp
          End Do
          ipC=0
          Do lAsh=1,naco
           il=0
           Do kash=1,naco
            Temp2(kash,lash)=0.0d0
            Do l=1,lbas*lcmp
             il=il+1
             Temp2(kash,lash)=Temp2(kash,lash)+
     &            Cl(ipc+l)*Temp1(il)
            End Do
           End Do
            ipC=ipC+lBas*lCmp
          End Do
*
          If (iShell(3).ne.ishell(4)) Then

          do iSPert=0,nIrrep-1
          If (pert(isPert)) Then
           rFact2=rFact*DBLE(iChtbl(ispert,nop(icnt)))
           llash=0
           k=abs(indgrd(icar,icnt,ispert))
           j=0
           Do lIrr=0,nIrrep-1
            Do lAsh=1,nAsh(lIrr)
             Do kIrr=0,lIrr
              irest=iEOR(iEOR(ioper(ispert),ioper(kirr)),ioper(lirr))
              kMax=nAsh(kIrr)
              If (kIrr.eq.lIrr) kMax=lAsh
              Do kAsh=1,kMax
               kk=kash+moip(kirr)
               ll=lash+moip(lirr)
               j=j+1
               Do jIrr=0,nIrrep-1
                rPj=DBLE(iChTbl(jIrr,nop(2)))
                 iirr=nropr(ieor(iOPER(jirr),irest))
                 rPij=rPj*DBLE(iChTbl(iIrr,nop(1)))*rfact2
                 buffer(ib,ic,jb,jc,iirr,j,k)=
     &               buffer(ib,ic,jb,jc,iirr,j,k)+
     &               rpij*Temp2(kk,ll)+
     &               rpij*Temp2(ll,kk)
               End Do
              End Do
             End Do
            End Do
           End Do
          End If
          End Do
          Else
*
          do iSPert=0,nIrrep-1
          If (pert(isPert)) Then
           rFact2=rFact*DBLE(iChtbl(ispert,nop(icnt)))
           llash=0
           k=abs(indgrd(icar,icnt,ispert))
           j=0
           Do lIrr=0,nIrrep-1
            Do lAsh=1,nAsh(lIrr)
             Do kIrr=0,lIrr
              irest=iEOR(iEOR(ioper(ispert),ioper(kirr)),ioper(lirr))
              kMax=nAsh(kIrr)
              If (kIrr.eq.lIrr) kMax=lAsh
              Do kAsh=1,kMax
               kk=kash+moip(kirr)
               ll=lash+moip(lirr)
               j=j+1
               Do jIrr=0,nIrrep-1
                rPj=DBLE(iChTbl(jIrr,nop(2)))
                iirr=nropr(ieor(iOPER(jirr),irest))
                rPij=rPj*DBLE(iChTbl(iIrr,nop(1)))*rfact2
                buffer(ib,ic,jb,jc,iirr,j,k)=
     &              buffer(ib,ic,jb,jc,iirr,j,k)+
     &              rpij*Temp2(kk,ll)
               End Do
              End Do
             End Do
            End Do
           End Do
          End If
          End Do
          End if
 1000     Continue
         End Do
        End Do
       End Do
      End Do
*
      kk=0
      Do kIrrep=0,nIrrep-1
       sfact=DBLE(ichtbl(kirrep,nop(3)))
       Do kAsh=1,nAsh(kIrrep)
        Do k=1,kcmp*kbas
         kk=kk+1
         Ck(kk)=Ck(kk)*sFact
        End Do
       End Do
      End Do
      kk=0
      Do kIrrep=0,nIrrep-1
       sfact=DBLE(ichtbl(kirrep,nop(4)))
       Do kAsh=1,nAsh(kIrrep)
        Do k=1,lcmp*lbas
         kk=kk+1
         Cl(kk)=Cl(kk)*sFact
        End Do
       End Do
      End Do

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nint)
         Call Unused_real(Temp3)
         Call Unused_real_array(MOInt)
         Call Unused_real(Ci)
         Call Unused_integer(nCi)
         Call Unused_real(Cj)
         Call Unused_integer(nCj)
         Call Unused_real_array(D)
         Call Unused_integer_array(iao)
         Call Unused_real(Tempi)
      End If
      End
