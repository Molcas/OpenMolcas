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
      SubRoutine r2elint_ns(rKappa,rMO1,rmo2,FockI,FockA,nF,
     &                   iDSym,sign,Fact,jspin)
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
      Logical lFI,lFA,lMo
      Parameter ( One = 1.0d0 )
      Real*8 rKappa(nDens2),rMO1(nMba),rmo2(*),FockI(nDens2),
     &       FockA(nDens2)
      Dimension rdum(1)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
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
      Call GetMem('Scr  ','ALLO','REAL',ipScr,imem)
      Call GetMem('Temp2','ALLO','REAL',ipTmp2,nDens22)
      Call GetMem('Temp3','ALLO','REAL',ipT3,nDens22)
      Call GetMem('Temp4','ALLO','REAL',ipT4,nDens22)
      Call GetMem('DIL  ','ALLO','REAL',ipDIL,nDens2)
      Call GetMem('DI   ','ALLO','REAL',ipDI,nCMO)
      Call GetMem('DIR  ','ALLO','REAL',ipDIR,nDens2)
      Call GetMem('FociI','ALLO','REAL',ipFI,ndens2)
      Call GetMem('FociI','ALLO','REAL',ipFI1,ndens2)
      Call GetMem('kappa','ALLO','REAL',ipK1,ndens2)

      call dcopy_(ndens2,[0.0d0],0,FockI,1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipFI),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipFI1),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipK1),1)
      call dcopy_(ndens2,[0.0d0],0,FockA,1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDir),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDil),1)
      call dcopy_(nCMO,[0.0d0],0,Work(ipDi),1)
      lFI=.true.
      lFa=.false.
      lMo=.false.
      If (iMethod.eq.2) Then
         Call GetMem('DAL  ','ALLO','REAL',ipDAL,nDens2)
         Call GetMem('DAR  ','ALLO','REAL',ipDAR,nDens2)
         Call GetMem('DA   ','ALLO','REAL',ipDA,nCMO)
         Call GetMem('FockA','ALLO','REAL',ipFA,nDens2)
         Call GetMem('FockA','ALLO','REAL',ipFA1,nDens2)
         lFa=.true.
         lMo=.true.
         call dcopy_(ndens2,[0.0d0],0,Work(ipFA),1)
         call dcopy_(ndens2,[0.0d0],0,Work(ipFA1),1)
         call dcopy_(ndens2,[0.0d0],0,Work(ipDar),1)
         call dcopy_(ndens2,[0.0d0],0,Work(ipDal),1)
         call dcopy_(nCMO,[0.0d0],0,Work(ipDa),1)
      Else
         ipDAL = ip_Dummy
         ipDAR = ip_Dummy
         ipDA  = ip_Dummy
         ipFA  = ip_Dummy
         ipFA1 = ip_Dummy
      End If
c THIS IS THE ENTIRE DENSITY FOR MULTICONF
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
          ip2=itri(iA,jA)
          Work(ipDA+ip-1)=Work(ipG1t+ip2-1)
         End Do
        End Do
       End Do
      End If

*
*     Construct {kappa,(ij|kl)} and the fock matrix contribution
*     from one index tranformation of contracted indexes
*
*      If (iDsym.eq.2) Then
*           Call RecPrt('Work(ipFI)',' ',Work(ipFI),ndens2,1)
*           Call RecPrt('Work(ipDI)',' ',Work(ipDI),ndens2,1)
*           Stop 10
*      End If
      FacR=Fact
      Call Read2_ns(rmo1,rmo2,
     &           FockI,FockA,
     &           Work(ipT1),imem,Work(ipScr),Work(ipTmp2),
     &           Work(ipT3),Work(ipT4),
     &           Work(ipDIR),Work(ipDIL),
     &           Work(ipDI),
     &           Work(ipDAR),Work(ipDAL),
     &           Work(ipDA),
     &           rkappa,idsym,Sign,Facr,jSpin,
     &           lFA,lfi,lMo,
     &           Work(ipCMO))
*      If (iDsym.eq.2) Then
*           Call RecPrt('Work(ipFI)',' ',Work(ipFI),ndens2,1)
*           Call RecPrt('Work(ipDI)',' ',Work(ipDI),ndens2,1)
*           Stop 10
*      End If
      Do iS=1,nsym
      js=ieor(idsym-1,is-1)+1
      If (nbas(is)*nbas(js).ne.0)
     &Call DGETMO(rkappa(ipmat(is,js)),nbas(is),nbas(is),nbas(js),
     &            Work(ipK1-1+ipmat(js,is)),nbas(js))
      End Do
      Call DSCAL_(ndens2,-1.0d0,Work(ipK1),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDir),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDil),1)
      If (imethod.eq.2) Then
      call dcopy_(ndens2,[0.0d0],0,Work(ipDar),1)
      call dcopy_(ndens2,[0.0d0],0,Work(ipDal),1)
      End If
      Call Read2_ns(rdum,rdum,
     &           Work(ipFi1),Work(ipFA1),
     &           Work(ipT1),imem,Work(ipScr),Work(ipTmp2),
     &           Work(ipT3),Work(ipT4),
     &           Work(ipDIR),Work(ipDIL),
     &           Work(ipDI),
     &           Work(ipDAR),Work(ipDAL),
     &           Work(ipDA),
     &           Work(ipK1),idsym,Sign,Facr,jSpin,
     &           lFA,lfi,.false.,
     &           Work(ipCMO))
*      If (iDsym.eq.2)
*     &      Call RecPrt('Work(ipFI1)',' ',Work(ipFI1),ndens2,1)
*
*      Calculate contribution from uncontracted indexes.
*
      Do iS=1,nSym
       jS=iEOr(iS-1,iDSym-1)+1
       If (nBas(iS)*nBas(jS).ne.0) Then
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            Work(ipFIMO+ipCM(iS)-1),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            Work(ipFIMO+ipCM(jS)-1),nBas(jS),
     &            One,FockI(ipMat(iS,jS)),nBas(is))
       If (iMethod.eq.2) Then
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            Work(ipFAMO+ipCM(iS)-1),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            Work(ipFAMO+ipCM(jS)-1),nBas(jS),
     &            One,FockA(ipMat(iS,jS)),nBas(is))
       End If
       End If
      End Do
*
      If (iMethod.eq.2) Then
         Call GetMem('FockA','FREE','REAL',ipFA,nF)
         Call GetMem('FockA','FREE','REAL',ipFA1,nF)
         Call GetMem('DAR  ','FREE','REAL',ipDAR,nDens2)
         Call GetMem('DA   ','FREE','REAL',ipDA,nDens2)
         Call GetMem('DAL  ','FREE','REAL',ipDAL,nDens2)
      End If
*
      Call GetMem('FociI','Free','REAL',ipFI,nF)
      Call GetMem('FociI','Free','REAL',ipFI1,nF)
      Call GETMEM('KAP1 ','FREE','REAL',ipK1,ndens2)
      Call GetMem('DIR  ','FREE','REAL',ipDIR,nDens2)
      Call GetMem('DIL  ','FREE','REAL',ipDIL,nDens22)
      Call GetMem('DI   ','FREE','REAL',ipDI,nDens2)
      Call GetMem('Temp4','FREE','REAL',ipT4,nDens22)
      Call GetMem('Temp3','FREE','REAL',ipT3,nDens22)
      Call GetMem('Temp2','FREE','REAL',ipTmp2,nDens22)
      Call GetMem('Scr  ','FREE','REAL',ipScr,imem)
      Call GetMem('Temp1','FREE','REAL',ipT1,imem)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
