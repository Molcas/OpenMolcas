!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************
      Subroutine DrvEMB_(nh1,KSDFT,Do_Grad,Grad,nGrad,DFTFOCK)
!***********************************************************************
!***********************************************************************
!** Orbital-Free Embedding calculation (gradients)                   ***
!**                                                                  ***
!**                                                                  ***
!** Author: F. Aquilante, Geneva Nov  2010                           ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
      use OFembed, only: Xsigma, dFMD
      Implicit Real*8 (a-h,o-z)
      External LSDA_emb, Checker
#include "real.fh"
#include "stdalloc.fh"
#include "debug.fh"
      Real*8 Grad(nGrad)
      Logical Do_Grad
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Character*16 NamRfil
      Real*8, Allocatable:: Grad_A(:), F_DFT(:,:), D_DS(:,:)
      Real*8, Allocatable:: Fcorr(:,:)
!
      Debug=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
      If (.not.Do_Grad) Then
         Call WarningMessage(2,'DrvEMB_: Do_Grad must be .true.')
         Call Abend()
      EndIf
      Call FZero(Grad,nGrad)
      Call mma_allocate(Grad_A,nGrad,Label='Grad_A')
      Grad_A(:)=Zero
!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys B (environment)             *
!                                                                      *
!***********************************************************************
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
!                                                                      *
!***********************************************************************
!                                                                      *
      nD=4
      Call mma_allocate(F_DFT,nh1,nD,Label='F_DFT')
      Call mma_allocate(D_DS,nh1,nD,Label='D_DS')
!
!---- Get the density matrix of the environment (rho_B)
!
      Call Get_iScalar('Multiplicity',kSpin)
      Call Get_D1ao(D_DS(1,1),nh1)
!     Call RecPrt('D_DS(1,1)',' ',D_DS(1,1),nh1,1)
!
!
!---- Get the spin density matrix of the environment
!
      If (kSpin.ne.1) Then
         Call Get_D1Sao(D_DS(1,2),nh1)
!        Call RecPrt('D1Sao',' ',D_DS(1,2),nh1,1)
      End If
!
!---- Compute alpha and beta density matrices of the environment
!
      nFckDim=2
      If (kSpin.eq.1) Then
         call dscal_(nh1,Half,D_DS(1,1),1)
         call dcopy_(nh1,D_DS(1,1),1,D_DS(1,2),1)
         nFckDim=1
      Else
         Do i = 1, nh1
            DTot=D_DS(i,1)
            DSpn=D_DS(i,2)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            D_DS(i,1)=d_Alpha
            D_DS(i,2)=d_Beta
         End Do
!      Call RecPrt('Da',' ',D_DS(1,1),nh1,1)
!      Call RecPrt('Db',' ',D_DS(1,2),nh1,1)
      End If
!

      If (KSDFT(1:4).eq.'NDSD') Then

         Call wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func_B,                    &
     &                   D_DS,nh1,nFckDim,                              &
     &                   Do_Grad,                                       &
     &                   Grad,nGrad,DFTFOCK)

         KSDFT(1:4)='LDTF' !set to Thomas-Fermi for subsequent calls
      EndIf
!                                                                      *
!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys A                           *
!                                                                      *
!***********************************************************************
      Call NameRun(NamRfil)    ! switch back RUNFILE name
!
!---- Get the density matrix for rho_A
!
      Call Get_D1ao(D_DS(1,3),nh1)
!     Call RecPrt('D_DS(1,3)',' ',D_DS(1,3),nh1,1)
!
      Call Get_iScalar('Multiplicity',iSpin)
      If (iSpin.eq.1 .and. kSpin.ne.1) Then
         Call WarningMessage(0,                                         &
     &     ' Non-singlet environment perturbation on singlet state!'//  &
     &     '  Spin-components of the OFE potential will be averaged. ' )
      EndIf
!
!---- Get the spin density matrix of A
!
      If (iSpin.ne.1) Then
         Call Get_D1Sao(D_DS(1,4),nh1)
!        Call RecPrt('D1Sao',' ',D_DS(1,4),nh1,1)
      End If
!
!---- Compute alpha and beta density matrices of subsystem A
!
      nFckDim=2
      If (iSpin.eq.1) Then
         call dscal_(nh1,Half,D_DS(1,3),1)
         call dcopy_(nh1,D_DS(1,3),1,D_DS(1,4),1)
         If (kSpin.eq.1) nFckDim=1
      Else
         Do i = 1, nh1
            DTot=D_DS(i,3)
            DSpn=D_DS(i,4)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            D_DS(i,3)=d_Alpha
            D_DS(i,4)=d_Beta
         End Do
!      Call RecPrt('Da',' ',D_DS(1,3),nh1,1)
!      Call RecPrt('Db',' ',D_DS(1,4),nh1,1)
      End If
!
      Call wrap_DrvNQ(KSDFT,F_DFT(1,3),nFckDim,Func_A,                  &
     &                D_DS(1,3),nh1,nFckDim,                            &
     &                Do_Grad,                                          &
     &                Grad_A,nGrad,DFTFOCK)

      Call daxpy_(nGrad,-1.0d0,Grad_A,1,Grad,1)
!
!  Fraction of correlation potential from A (cases: HF or Trunc. CI)
      If (dFMD.gt.0.0d0) Then
!
         Grad_A(:)=Zero
         Call mma_allocate(Fcorr,nh1,nFckDim,Label='Fcorr')

         Call cwrap_DrvNQ(KSDFT,F_DFT(1,3),nFckDim,Func_A,              &
     &                    D_DS(1,3),nh1,nFckDim,                        &
     &                    Do_Grad,                                      &
     &                    Grad_A,nGrad,DFTFOCK,Fcorr)

         Call get_dScalar('NAD dft energy',Energy_NAD)
         Fakt_ = Xlambda(abs(Energy_NAD),Xsigma)
         Call daxpy_(nGrad,Fakt_,Grad_A,1,Grad,1)

         Call mma_deallocate(Fcorr)
      End If
!
      Call mma_deallocate(Grad_A)
!
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
!
      Call wrap_DrvNQ('NUCATT_EMB',F_DFT,nFckDim,Func_X,                &
     &                D_DS(1,3),nh1,nFckDim,                            &
     &                Do_Grad,                                          &
     &                Grad,nGrad,DFTFOCK)
!
      Call NameRun(NamRfil)   ! switch back RUNFILE name
!
!***********************************************************************
!                                                                      *
!     Calculation on the supermolecule                                 *
!                                                                      *
!***********************************************************************
      nFckDim=2
      If (iSpin.eq.1 .and. kSpin.eq.1) Then
         nFckDim=1
         Call daxpy_(nh1,One,D_DS(1,3),1,D_DS(1,1),1)
      Else
         Call daxpy_(nh1,One,D_DS(1,3),1,D_DS(1,1),1)
         Call daxpy_(nh1,One,D_DS(1,4),1,D_DS(1,2),1)
      EndIf

      Call wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func_AB,                      &
     &                D_DS,nh1,nFckDim,                                 &
     &                Do_Grad,                                          &
     &                Grad,nGrad,DFTFOCK)
!
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(D_DS)
!
      Return
      End
