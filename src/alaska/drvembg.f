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
* Copyright (C) 2010, Francesco Aquilante                              *
************************************************************************
      SubRoutine DrvEMBg(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
* Object: driver for computation of gradient with respect to the OFE   *
*         energy (DFT-type contribution only)                          *
*                                                                      *
*                                                                      *
*         int{ [Vxc^nad(r) + Ts^nad(r) + Vnuc^B(r)] (drhoA/dRk) dr }   *
*                                                                      *
*                                                                      *
* Note:  rho=rho_A+rho_B  where the dmat elements of rho_B on the bsfs *
*                         of A have been annihilated. This also gives  *
*                         drhoA/dRk=drho/dRk (for Rk in A)             *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
*                                                                      *
*     Author: F. Aquilante, Dept. of Phys. Chem.                       *
*             University of Geneva, Switzerland                        *
************************************************************************
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
      Character*16  OFE_KSDFT
      COMMON  / OFembed_C / OFE_KSDFT
      Character Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
      Logical Do_Grad, King
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      Call qEnter('DrvEMBg')
      LuWr=6
*
      Call StatusLine(' Alaska:',' Computing OFembedding gradients')
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*     DFT-OFE    gradient                                              *
************************************************************************
*                                                                      *
      Do_Grad=.True.
      Call DrvEMB_(nDens,OFE_KSDFT,Do_Grad,Temp,nGrad,'SCF ')
*
      iEnd=1
 99   Continue
      If (OFE_KSDFT(iEnd:iEnd).eq.' ') Then
         iEnd = iEnd - 1
      Else
         iEnd = iEnd + 1
         Go To 99
      End If
      Label='DFT-OFE('//OFE_KSDFT(1:iEnd)//') contribution'
      jPrint=nPrint(112)
      If (jPrint.ge.15) Call PrGrad(Label,Temp,nGrad,ChDisp,5)
      If (king()) Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      If (iPrint.lt.6) Go To 777
      Write (LuWr,*)
      If (Grid_Type.eq.Moving_Grid) Then
         Write(LuWr,*)'DFT-OFE contribution computed for a moving grid.'
      Else
         Write(LuWr,*)'DFT-OFE contribution computed for a fixed grid.'
      End If
      Write (LuWr,*)
 777  Continue
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call CWTime(TCpu2,TWall2)
      Call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
      Call qExit('DrvEMBg')
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine DrvEMB_(nh1,KSDFT,Do_Grad,Grad,nGrad,DFTFOCK)
************************************************************************
************************************************************************
*** Orbital-Free Embedding calculation (gradients)                   ***
***                                                                  ***
***                                                                  ***
*** Author: F. Aquilante, Geneva Nov  2010                           ***
***                                                                  ***
************************************************************************
************************************************************************
      Implicit Real*8 (a-h,o-z)
      External LSDA_emb, Checker
#include "real.fh"
#include "WrkSpc.fh"
#include "debug.fh"
      Real*8 Grad(nGrad)
      Logical Do_Grad
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Character*16 NamRfil
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_R1/ Xsigma
      COMMON  / OFembed_R2/ dFMD
*
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('DrvEMB_')
      If (.not.Do_Grad) Then
         Call WarningMessage(2,'DrvEMB_: Do_Grad must be .true.')
         Call Abend()
      EndIf
      Call FZero(Grad,nGrad)
      Call GetMem('Grad_A','Allo','Real',ip_Grad_A,nGrad)
      Call FZero(Work(ip_Grad_A),nGrad)
************************************************************************
*                                                                      *
*     Setup of density matrices for subsys B (environment)             *
*                                                                      *
************************************************************************
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
*                                                                      *
************************************************************************
*                                                                      *
      nD=4
      lFck=nh1*nD
      Call Allocate_Work(ipF_DFT,lFck)
      ipFA_DFT=ipF_DFT+2*nh1
      l_D_DS=nh1*nD
      Call GetMem('D-DS','Allo','Real',ip_D_DS,l_D_DS)
      ipA_D_DS=ip_D_DS+2*nh1
*
*---- Get the density matrix of the environment (rho_B)
*
      Call Get_iScalar('Multiplicity',kSpin)
      Call Get_D1ao(ipD1ao,nDens)
      If (nDens.ne.nh1) Then
         Call WarningMessage(2,'DrvEMB_: nDens.ne.nh1')
         Write (6,*) 'nDens=',nDens
         Write (6,*) 'nh1  =',nh1
         Call Abend()
      End If
      call dcopy_(nh1,Work(ipD1ao),1,Work(ip_D_DS),1)
*     Call RecPrt('D1ao',' ',Work(ipD1ao),nh1,1)
*
      Call GetMem('Dens','Free','Real',ipD1ao,nDens)
*
*---- Get the spin density matrix of the environment
*
      If (kSpin.ne.1) Then
         Call Get_D1Sao(ipD1Sao,nDens)
*        Call RecPrt('D1Sao',' ',Work(ipD1Sao),nh1,1)
         call dcopy_(nh1,Work(ipD1Sao),1,Work(ip_D_DS+nh1),1)
         Call GetMem('Dens','Free','Real',ipD1Sao,nDens)
      End If
*
*---- Compute alpha and beta density matrices of the environment
*
      nFckDim=2
      If (kSpin.eq.1) Then
         call dscal_(nh1,Half,Work(ip_D_DS),1)
         call dcopy_(nh1,Work(ip_D_DS),1,Work(ip_D_DS+nh1),1)
         nFckDim=1
      Else
         Do i = 1, nh1
            DTot=Work(ip_D_DS+i-1)
            DSpn=Work(ip_D_DS+i-1+nh1)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            Work(ip_D_DS+i-1)=    d_Alpha
            Work(ip_D_DS+i-1+nh1)=d_Beta
         End Do
*      Call RecPrt('Da',' ',Work(ip_D_DS),nh1,1)
*      Call RecPrt('Db',' ',Work(ip_D_DS+nh1),nh1,1)
      End If
*

      If (KSDFT(1:4).eq.'NDSD') Then

         Call wrap_DrvNQ(KSDFT,Work(ipF_DFT),nFckDim,Func_B,
     &                   Work(ip_D_DS),nh1,nFckDim,
     &                   Do_Grad,
     &                   Grad,nGrad,DFTFOCK)

         KSDFT(1:4)='LDTF' !set to Thomas-Fermi for subsequent calls
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*     Setup of density matrices for subsys A                           *
*                                                                      *
************************************************************************
      Call NameRun(NamRfil)    ! switch back RUNFILE name
*
*---- Get the density matrix for rho_A
*
      Call Get_D1ao(ipD1ao,nDens)
      If (nDens.ne.nh1) Then
         Call WarningMessage(2,'DrvEMB_: nDens.ne.nh1')
         Write (6,*) 'nDens=',nDens
         Write (6,*) 'nh1  =',nh1
         Call Abend()
      End If
      call dcopy_(nh1,Work(ipD1ao),1,Work(ipA_D_DS),1)
*     Call RecPrt('D1ao',' ',Work(ipD1ao),nh1,1)
*
      Call GetMem('Dens','Free','Real',ipD1ao,nDens)
*
      Call Get_iScalar('Multiplicity',iSpin)
      If (iSpin.eq.1 .and. kSpin.ne.1) Then
         Call WarningMessage(0,
     &     ' Non-singlet environment perturbation on singlet state!'//
     &     '  Spin-components of the OFE potential will be averaged. ' )
      EndIf
*
*---- Get the spin density matrix of A
*
      If (iSpin.ne.1) Then
         Call Get_D1Sao(ipD1Sao,nDens)
*        Call RecPrt('D1Sao',' ',Work(ipD1Sao),nh1,1)
         call dcopy_(nh1,Work(ipD1Sao),1,Work(ipA_D_DS+nh1),1)
         Call GetMem('Dens','Free','Real',ipD1Sao,nDens)
      End If
*
*---- Compute alpha and beta density matrices of subsystem A
*
      nFckDim=2
      If (iSpin.eq.1) Then
         call dscal_(nh1,Half,Work(ipA_D_DS),1)
         call dcopy_(nh1,Work(ipA_D_DS),1,Work(ipA_D_DS+nh1),1)
         If (kSpin.eq.1) nFckDim=1
      Else
         Do i = 1, nh1
            DTot=Work(ipA_D_DS+i-1)
            DSpn=Work(ipA_D_DS+i-1+nh1)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            Work(ipA_D_DS+i-1)=    d_Alpha
            Work(ipA_D_DS+i-1+nh1)=d_Beta
         End Do
*      Call RecPrt('Da',' ',Work(ipA_D_DS),nh1,1)
*      Call RecPrt('Db',' ',Work(ipA_D_DS+nh1),nh1,1)
      End If
*
      Call wrap_DrvNQ(KSDFT,Work(ipFA_DFT),nFckDim,Func_A,
     &                Work(ipA_D_DS),nh1,nFckDim,
     &                Do_Grad,
     &                Work(ip_Grad_A),nGrad,DFTFOCK)

      Call daxpy_(nGrad,-1.0d0,Work(ip_Grad_A),1,Grad,1)
*
*  Fraction of correlation potential from A (cases: HF or Trunc. CI)
      If (dFMD.gt.0.0d0) Then
*
         Call FZero(Work(ip_Grad_A),nGrad)
         Call GetMem('Fcorr','Allo','Real',ipFc,nh1*nFckDim) !dummy

         Call cwrap_DrvNQ(KSDFT,Work(ipFA_DFT),nFckDim,Func_A,
     &                    Work(ipA_D_DS),nh1,nFckDim,
     &                    Do_Grad,
     &                    Work(ip_Grad_A),nGrad,DFTFOCK,Work(ipFc))

         Call get_dScalar('NAD dft energy',Energy_NAD)
         Fakt_ = Xlambda(abs(Energy_NAD),Xsigma)
         Call daxpy_(nGrad,Fakt_,Work(ip_Grad_A),1,Grad,1)

         Call GetMem('Fcorr','Free','Real',ipFc,nh1*nFckDim)
      End If
*
      Call GetMem('Grad_A','Free','Real',ip_Grad_A,nGrad)
*
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
*
      Call wrap_DrvNQ('NUCATT_EMB',Work(ipF_DFT),nFckDim,Func_X,
     &                Work(ipA_D_DS),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)
*
      Call NameRun(NamRfil)   ! switch back RUNFILE name
*
************************************************************************
*                                                                      *
*     Calculation on the supermolecule                                 *
*                                                                      *
************************************************************************
      nFckDim=2
      If (iSpin.eq.1 .and. kSpin.eq.1) Then
         nFckDim=1
         Call daxpy_(nh1,One,Work(ipA_D_DS),1,Work(ip_D_DS),1)
      Else
         Call daxpy_(nh1,One,Work(ipA_D_DS),1,Work(ip_D_DS),1)
         Call daxpy_(nh1,One,Work(ipA_D_DS+nh1),1,Work(ip_D_DS+nh1),1)
      EndIf

      Call wrap_DrvNQ(KSDFT,Work(ipF_DFT),nFckDim,Func_AB,
     &                Work(ip_D_DS),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)
*
      Call Free_Work(ipF_DFT)
      Call GetMem('D-DS','Free','Real',ip_D_DS,l_D_DS)
      Call QExit('DrvEMB_')
*
      Return
      End
