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
      SubRoutine r2elint_sp(rKappa,rMO1,rmo2,FockI,FockA,nF,
     &                   iDSym,sign,Fact,jspin,D,FA)
*
************************************************************************
*
*     Constructs the one index transformed Fock-matrixes
*     and (pj|kl).
*     rKappa : the transformation matrix
*     iDSym  : Symmetry of transformation
*     sign   : 1:real -1:complex
*     jspin  : 0:singlet 1:triplet
*
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "spin.fh"
      Logical lFI,lFA,lMo
      Parameter ( One = 1.0d0 )
      Real*8 rKappa(nDens2),rMO1(nMba),rmo2(*),FockI(nDens2),
     &       FockA(nDens2),D(*),FA(*)
*
      irec(i,j)=i+nna*(j-1)
      ndens22=ndens2
      iAM=0
      iBM=0
      Do i=1,nSym
       Do j=1,nSym
        nDens22=Max(nDens22,nBas(i)*nBas(j))
       End do
       iAM=Max(iAM,nAsh(i)+nIsh(i))
       iBM=Max(iBM,nBAs(i))
      End Do
      imem=nDens22
      Call GetMem('Temp1','ALLO','REAL',ipT1,imem)
      Call GetMem('Temp2','ALLO','REAL',ipTmp2,nDens22)
      Call GetMem('Temp3','ALLO','REAL',ipT3,nDens22)
      Call GetMem('Temp4','ALLO','REAL',ipT4,nDens22)
      Call GetMem('DIL  ','ALLO','REAL',ipDIL,nDens2)
      Call GetMem('DI   ','ALLO','REAL',ipDI,nCMO)
      Call GetMem('DIR  ','ALLO','REAL',ipDIR,nDens2)
      Call GetMem('FociI','ALLO','REAL',ipFI,ndens2)
      call dcopy_(ndens2,[0.0d0],0,FockI,1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipFI),1)
C     Call GetMem('FockI','CHECK','REAL',ipFI,ndens2)
      call dcopy_(ndens2,[0.0d0],0,FockA,1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDir),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDil),1)
C     Call GetMem('FockI','CHECK','REAL',ipFI,ndens2)
      lFI=.true.
      lFa=.false.
      lMo=.false.
      If (iMethod.eq.2) Then
         Call GetMem('DAL  ','ALLO','REAL',ipDAL,nDens2)
         Call GetMem('DAR  ','ALLO','REAL',ipDAR,nDens2)
         Call GetMem('DA   ','ALLO','REAL',ipDA,nCMO)
         Call GetMem('FockA','ALLO','REAL',ipFA,nDens2)
         lFa=.true.
         lMo=.true.
      Else
         ipDAL = ip_Dummy
         ipDAR = ip_Dummy
         ipDA  = ip_Dummy
         ipFA  = ip_Dummy
      End If
      Do iS=1,nSym
      Do iB=1,nIsh(iS)
      ip=ipCM(iS)+(ib-1)*nBas(is)+ib-1
      Work(ipDI+ip-1)=2.0d0
      End Do
      End Do
      If (iMethod.eq.2) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
          iA=nA(is)+ib
          jA=nA(is)+jb
          ip2=irec(iA,jA)
          Work(ipDA+ip-1)=D(ip2)
         End Do
        End Do
       End Do
      End If

*     iiii=ipDA
*     ipDA=ipDI
*     ipFAMO=ipFIMO
*
*     Construct {kappa,(ij|kl)} and the fock matrix contribution
*     from one index tranformation of contracted indexes
*
      FacR=Fact
      Call Read2_2(rmo1,rmo2,
     &           Work(ipFI),Work(ipFA),
     &           Work(ipT1),imem,Work(ipTmp2),
     &           Work(ipT3),Work(ipT4),
     &           Work(ipDIR),Work(ipDIL),
     &           Work(ipDI),
     &           Work(ipDAR),Work(ipDAL),
     &           Work(ipDA),
     &           rkappa,idsym,Sign,Facr,jSpin,
     &           lFA,lfi,lMo,
     &           Work(ipCMO))
*
*      Calculate contribution from uncontracted indexes.
*
      Do iS=1,nSym
       jS=iEOr(iS-1,iDSym-1)+1
       If (nBas(iS)*nBas(jS).ne.0) Then
       Call DGEMM_('T','N',
     &             nBas(iS),nBas(jS),nBas(iS),
     &             1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &             Work(ipFI-1+ipMat(iS,jS)),nBas(iS),
     &             0.0d0,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            Work(ipFIMO+ipCM(iS)-1),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            Work(ipFIMO+ipCM(jS)-1),nBas(jS),
     &            One,FockI(ipMat(iS,jS)),nBas(is))
       If (iMethod.eq.2) Then
         Call DGEMM_('T','N',
     &               nBas(iS),nBas(jS),nBas(iS),
     &               1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &               Work(ipFA-1+ipMat(iS,jS)),nBas(iS),
     &               0.0d0,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            FA(ipCM(iS)),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            FA(ipCM(jS)),nBas(jS),
     &            One,FockA(ipMat(iS,jS)),nBas(is))
       End If
       End If
      End Do
      If (iMethod.eq.iCASSCF) Then
*       Call Dyax(nmba,1.0d0,Work(ipMT1),1,rmo,1)
*       Call DaXpY_(nmba,1.0d0,Work(ipMT2),1,rmo,1)
        Call GetMem('FockA','FREE','REAL',ipFA,nF)
        Call GetMem('DAR  ','FREE','REAL',ipDAR,nDens2)
        Call GetMem('DA   ','FREE','REAL',ipDA,nDens2)
        Call GetMem('DAL  ','FREE','REAL',ipDAL,nDens2)
      End If
      Call GetMem('FociI','Free','REAL',ipFI,nF)
      Call GetMem('DIR  ','FREE','REAL',ipDIR,nDens2)
      Call GetMem('DIL  ','FREE','REAL',ipDIL,nDens22)
      Call GetMem('DI   ','FREE','REAL',ipDI,nDens2)
      Call GetMem('Temp4','FREE','REAL',ipT4,nDens22)
      Call GetMem('Temp3','FREE','REAL',ipT3,nDens22)
      Call GetMem('Temp2','FREE','REAL',ipTmp2,nDens22)
      Call GetMem('Temp1','FREE','REAL',ipT1,imem)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
