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
      Use Basis_Info, only: dbsc, nCnttp, nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External OvrGrd, KneGrd, NAGrd, PrjGrd, M1Grd, M2Grd, SROGrd,
     &         WelGrd, XFdGrd, RFGrd, PCMGrd, PPGrd, COSGrd, FragPGrd
      External OvrMmG, KneMmG, NAMmG, PrjMmG, M1MmG, M2MmG, SROMmG,
     &         WelMmg, XFdMmg, RFMmg, PCMMmg, PPMmG, FragPMmG
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"
      Integer, Allocatable:: lOper(:)
      Real*8, Allocatable:: Coor(:,:), D_Var(:)

      Character Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
CAOM<
      parameter (lforce=20)
      real*8 force(lforce)
      common /finfld/force
      External MltGrd,MltMmG
CAOM>
      Logical DiffOp, lECP, lPP, lFAIEMP
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
      Call StatusLine(' Alaska:',' Computing 1-el OFE gradients')
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
      lECP=.False.
      lPP =.False.
      lFAIEMP=.False.
      Do i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
         lPP  = lPP  .or. dbsc(i)%nPP.ne.0
         lFAIEMP = LFAIEMP .or. dbsc(i)%Frag
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

      Call mma_allocate(D_Var,nDens,Label='D_Var')
      Call Get_D1ao_Var(D_var,nDens)
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=1
         Do iIrrep = 0, nIrrep - 1
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',D_Var(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
* Annihilate all the components of rho_B in the bsfs of the A subsystem
*
*SVC: fixed according to instructions from Francesco,
* as embedding should not deal with symmetry
      Call Annihil_rho(D_var,nBas(0))
*
      Call NameRun(NamRfil)   ! switch RUNFILE name
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     lOper(:): lOper of each component of the operator
*
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call mma_allocate(Coor,3,nComp)
      Call mma_allocate(lOper,nComp,Label='lOper')
      Coor(:,:)=Zero
      lOper(:) = 1
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
      Call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Coor,
     &             D_Var,nDens,lOper,nComp,nOrdOp,Label)

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
         Call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M1 Operator contribution'
         Call OneEl_g( M1Grd, M1MmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M2 Operator contribution'
         Call OneEl_g( M2Grd, M2MmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The SR Operator contribution'
         Call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
      If (lPP) Then
         Label = ' The Pseudo Potential contribution'
         Call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Coor,
     &                 D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
*                                                                      *
************************************************************************
*
      If(lFAIEMP) Then
        DiffOp = .True.
        Label = ' The FAIEMP Projection Operator Contribution'
        Call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Coor,
     &               D_Var,nDens,lOper,nComp,nOrdOp,Label)
        Call DaXpY_(nGrad,One,Temp,1,Grad,1)
        Call DrvG_FAIEMP(Grad,Temp,nGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(lOper)
      Call mma_deallocate(Coor)
*                                                                      *
************************************************************************
*                                                                      *
*...  Epilogue, end
*
      Call mma_deallocate(D_Var)
*
      Call Free_iSD()
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Annihil_rho(Dmat,nBas)

      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 Dmat(*)
      Integer nBas
      Character*(LENIN4) Name(mxBas)
      Integer, Allocatable:: nBas_per_Atom(:), nBas_Start(:)
      Real*8, Allocatable:: Charge_B(:)
*
************************************************************************
*
      Call Get_iScalar('Unique atoms',nAtoms)

      If (nAtoms.lt.1) Then
         Write(6,'(A,I9)') 'nUniqAt =',nAtoms
         Call Abend()
      End If
*
      Call mma_allocate(nBas_per_Atom,nAtoms,Label='nBpA')
      Call mma_allocate(nBas_Start,nAtoms,Label='nB_Start')

      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas)

      Call BasFun_Atom(nBas_per_Atom,nBas_Start,
     &                 Name,nBas,nAtoms,.false.)
*
      Call mma_allocate(Charge_B,nAtoms,Label='Charge_B')
      Call Get_dArray('Nuclear charge',Charge_B,nAtoms)
*
      ZA=0.0d0
      iAt=1
      Do while (iAt.le.nAtoms .and. ZA.eq.0.0d0)
         ZA = Charge_B(iAt)
         iAt = iAt + 1
      End Do
      iAt_B=iAt-1  ! start of atoms of subsystem B
      Call mma_deallocate(Charge_B)
*
      If (iAt_B.eq.1) Then ! subsystem B comes first
         ZB=1.0d0
         nAt_B=1
         Do while (nAt_B.le.nAtoms .and. ZB.gt.0.0d0)
            ZB = Charge_B(nAt_B)
            nAt_B = nAt_B + 1
         End Do
         nAt_B=nAt_B-1  ! end of atoms of subsystem B
         nBas_B = nBas_Start(nAt_B) - 1
         Do j=nBas_B,nBas-1
            jj=j*(j+1)/2
            Do i=1,j
               ijj=i+jj
               Dmat(ijj)=0.0d0
            End Do
         End Do
      Else
         nBas_A = nBas_Start(iAt_B) - 1
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
      Call mma_deallocate(nBas_Start)
      Call mma_deallocate(nBas_per_Atom)
*
*  Annihilated density written to runfile for use in Coulomb gradients
*
      Length=nBas*(nBas+1)/2
      Call Put_D1ao_Var(Dmat,Length)

      Return
      End
