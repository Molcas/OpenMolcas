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
      SubRoutine Drvh1_EMB(Grad,Temp,nGrad)
      Use Basis_Info, only: dbsc, nCnttp
      Implicit Real*8 (A-H,O-Z)
      External OvrGrd, KneGrd, NAGrd, PrjGrd, M1Grd, M2Grd, SROGrd,
     &         WelGrd, XFdGrd, RFGrd, PCMGrd, PPGrd, COSGrd, FragPGrd
      External OvrMmG, KneMmG, NAMmG, PrjMmG, M1MmG, M2MmG, SROMmG,
     &         WelMmg, XFdMmg, RFMmg, PCMMmg, PPMmG, FragPMmG
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"
      Character Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
CAOM<
      parameter (lforce=20)
      real*8 force(lforce)
      common /finfld/force
      External MltGrd,MltMmG
CAOM>
      Logical DiffOp, lECP, lPP
      Character*16 NamRfil
*
*-----Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
*...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      Call CWTime(TCpu1,TWall1)
      Call qEnter('Drvh1_emb')
      Call StatusLine(' Alaska:',' Computing 1-el OFE gradients')
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
      lECP=.False.
      lPP =.False.
      Do i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
         lPP  = lPP  .or. dbsc(i)%nPP.ne.0
      End Do
*
*---- Allocate memory for density matrices
*
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis

      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name

      Call Get_D1ao_Var(ipD_var,Length)
      If ( length.ne.nDens ) Then
         Call WarningMessage(2,'Error in Drvh1_emb')
         Write (6,*) 'Drvh1_emb: length.ne.nDens'
         Write (6,*) 'length,nDens=',length,nDens
         Call Abend()
      End If
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=ipD_Var
         Do iIrrep = 0, nIrrep - 1
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',Work(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
* Annihilate all the components of rho_B in the bsfs of the A subsystem
*
*SVC: fixed according to instructions from Francesco,
* as embedding should not deal with symmetry
      Call Annihil_rho(Work(ipD_var),nBas(0))
*
      Call NameRun(NamRfil)   ! switch RUNFILE name
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     Work(ip1): lOper of each component of the operator
*
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call GetMem('Coor','Allo','Real',ipC,3*nComp)
      Call GetMem('lOper','Allo','Inte',ip1,nComp)
      call dcopy_(nComp*3,[Zero],0,Work(ipC),1)
      iWork(ip1) = 1
*
************************************************************************
*3)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the nuclear attraction integrals.                    *
*                                                                      *
************************************************************************
*
      DiffOp = .True.
      Label = ' The Nuclear Attraction Contribution'
      Call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Work(ipC),
     &           Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)

      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
************************************************************************
*4)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the ECP integrals.                                   *
*                                                                      *
************************************************************************
*
      If (lECP) Then
         DiffOp = .True.
         Label = ' The Projection Operator contribution'
         Call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M1 Operator contribution'
         Call OneEl_g( M1Grd, M1MmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M2 Operator contribution'
         Call OneEl_g( M2Grd, M2MmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The SR Operator contribution'
         Call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
      If (lPP) Then
         Label = ' The Pseudo Potential contribution'
         Call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
*                                                                      *
************************************************************************
*
      If(lFAIEMP) Then
        DiffOp = .True.
        Label = ' The FAIEMP Projection Operator Contribution'
        Call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Work(ipC),
     &               Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
        Call DaXpY_(nGrad,One,Temp,1,Grad,1)
        Call DrvG_FAIEMP(Grad,Temp,nGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('lOper','Free','Inte',ip1,nComp)
      Call GetMem('Coor','Free','Real',ipC,3*nComp)
*                                                                      *
************************************************************************
*                                                                      *
*...  Epilogue, end
*
      Call GetMem('D0  ','Free','Real',ipD_Var,nDens)
*
      Call Free_iSD()
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Call qExit('Drvh1_emb')
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Annihil_rho(Dmat,nBas)

      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      Real*8 Dmat(*)
      Integer nBas
      Character*(LENIN4) Name(mxBas)
*
************************************************************************
*     Statement function for Charges of subsystem B
      Charge_B(i) = Work(ip_ChargeB+i-1)
************************************************************************
*
      Call Get_iScalar('Unique atoms',nAtoms)

      If (nAtoms.lt.1) Then
         Write(6,'(A,I9)') 'nUniqAt =',nAtoms
         Call Abend()
      End If
*
      Call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,nAtoms)
      Call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,nAtoms)

      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas)

      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 Name,nBas,nAtoms,.false.)
*
      Call GetMem('ChargeB','Allo','Real',ip_ChargeB,nAtoms)
      Call Get_dArray('Nuclear charge',Work(ip_ChargeB),nAtoms)
*
      ZA=0.0d0
      iAt=1
      Do while (iAt.le.nAtoms .and. ZA.eq.0.0d0)
         ZA = Charge_B(iAt)
         iAt = iAt + 1
      End Do
      iAt_B=iAt-1  ! start of atoms of subsystem B
      Call GetMem('ChargeB','Free','Real',ip_ChargeB,nAtoms)
*
      If (iAt_B.eq.1) Then ! subsystem B comes first
         ZB=1.0d0
         nAt_B=1
         Do while (nAt_B.le.nAtoms .and. ZB.gt.0.0d0)
            ZB = Charge_B(nAt_B)
            nAt_B = nAt_B + 1
         End Do
         nAt_B=nAt_B-1  ! end of atoms of subsystem B
         nBas_B = iWork(ip_nBas_Start+nAt_B-1) - 1
         Do j=nBas_B,nBas-1
            jj=j*(j+1)/2
            Do i=1,j
               ijj=i+jj
               Dmat(ijj)=0.0d0
            End Do
         End Do
      Else
         nBas_A = iWork(ip_nBas_Start+iAt_B-1) - 1
         nAA=nBas_A*(nBas_A+1)/2
         Call FZero(Dmat,nAA)
         Do j=nBas_A,nBas-1
            jj=j*(j+1)/2
            Do i=1,nBas_A
               ijj=i+jj
               Dmat(ijj)=0.0d0
            End Do
         End Do
      EndIf
*
      Call GetMem('nB_Start','Free','Inte',ip_nBas_Start,nAtoms)
      Call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,nAtoms)
*
*  Annihilated density written to runfile for use in Coulomb gradients
*
      Length=nBas*(nBas+1)/2
      Call Put_D1ao_Var(Dmat,Length)
      Return
      End
