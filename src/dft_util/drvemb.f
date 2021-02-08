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
* Copyright (C) 2010,2012,2017, Francesco Aquilante                    *
*               2015,2017, Alexander Zech                              *
************************************************************************
      Subroutine DrvEMB(h1,D,RepNuc,nh1,
     &                  KSDFT,ExFac,Do_Grad,Grad,nGrad,
     &                  D1I,D1A,nD1,DFTFOCK)
************************************************************************
************************************************************************
*** Orbital-Free Embedding calculation                               ***
***                                                                  ***
*** Method:                                                          ***
***     T. A. Wesolowski, A. Warshel, J. Phys. Chem. 97 (1993) 8050. ***
***                                                                  ***
*** NDSD potential:                                                  ***
***     J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski,       ***
***                               J. Chem. Phys.  129 (2008) 074107. ***
***                                                                  ***
*** Embedding multi-determinantal wfs:                               ***
***     T. A. Wesolowski, Phys. Rev.A. 77 (2008) 012504.             ***
***                                                                  ***
***                                                                  ***
*** Embedding Hartree-Fock wf:                                       ***
***     F. Aquilante, T. A. Wesolowski                               ***
***                       J. Chem. Phys. 135 (2011) 084120.          ***
***                                                                  ***
***                                                                  ***
*** Author: F. Aquilante, Geneva July 2010                           ***
***                                                                  ***
***                       (last update: Feb 2012)                    ***
***                                                                  ***
************************************************************************
************************************************************************
      use OFembed, only: OFE_first, Xsigma, dFMD
      use OFembed, only: Func_AB,Func_A,Func_B,Energy_NAD,
     &                   V_Nuc_AB,V_Nuc_BA,V_emb
      Implicit Real*8 (a-h,o-z)
      External LSDA_emb, Checker
#include "real.fh"
#include "stdalloc.fh"
#include "debug.fh"
      Real*8 h1(nh1), D(nh1,2), Grad(nGrad)
      Real*8 D1I(nD1),D1A(nD1)
      Logical Do_Grad
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Character*16 NamRfil
      Real*8 Vxc_ref(2)
*
*     Real*8 Func_A_TF, Func_B_TF, Func_AB_TF, TF_NAD
      Real*8 Func_A_TF, Func_B_TF
      Logical is_rhoA_on_file
      Real*8 Xlambda
      External Xlambda
*
      Real*8, Allocatable:: D_DS(:,:), F_DFT(:,:), Fcorr(:,:), TmpA(:)
#ifdef _NOT_USED_
      Real*8, Allocatable:: Vemb(:), D1ao_x(:)
#endif
*
      Debug=.False.
      is_rhoA_on_file = .False.
*                                                                      *
************************************************************************
*                                                                      *
      Call Setup_iSD()
      If (Do_Grad) Call FZero(Grad,nGrad)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _NOT_USED_
*     --- Section to calculate Nonelectr. V_emb with current density
*     Temporarily turned off (clean output)
       If (.not.OFE_first) then
          Call mma_allocate(D1ao_y,nh1)
          Call Get_NameRun(NamRfil) ! save the old RUNFILE name
          Call NameRun('AUXRFIL')   ! switch RUNFILE name
          Call mma_allocate(Vemb,nh1,label='Vemb')
          Call Get_dArray('dExcdRa', Vemb, nh1)
          Call mma_allocate(TmpA,nh1,Label='TmpA')
          Call Get_dArray('Nuc Potential',TmpA,nh1)
*    Substract V_nuc_B
          Call daxpy_(nh1,-One,TmpA,1,Vemb,1)
*    Calculate nonelectr. V_emb with current Density
          Ynorm=dDot_(nh1,WD1ao_y,1,D1ao_y,1)
          V_emb_x=dDot_(nh1,Vemb,1,D1ao_y,1)
          Write (6,'(A,F19.10,4X,A,F10.5)')
     &          'Nonelectr. Vemb w. current density: ', V_emb_x,
     &          'Y_Norm = ', Ynorm
          Call mma_deallocate(D1ao_y)
*    Get rho_A_ref
          Call NameRun('PRERFIL')
          Call mma_allocate(D1ao_x,nDens,Label='D1ao_x')
          Call get_dArray('D1ao',ipD1ao_x,nDens)
          Xnorm=dDot_(nh1,D1ao_x,1,D1ao_x,1)
          V_emb_x_ref=dDot_(nh1,Vemb,1,pD1ao_x,1)
          Write (6,'(A,F19.10,4X,A,F10.5)')
     &          'Nonelectr. Vemb w.    ref. density: ', V_emb_x_ref,
     &          'X_Norm = ', Xnorm
          Call VEMB_Exc_states(Vemb,nh1,KSDFT,Func_B)
          Call mma_deallocate(TmpA)
          Call mma_deallocate(D1ao_x)
          Call mma_dealloacte(Vemb)
          Call NameRun(NamRfil)     ! switch back to RUNFILE
      End If
*     --- Section End
#endif
      Call f_Inquire('PRERFIL',is_rhoA_on_file) ! rho_A from file
      If (is_rhoA_on_file .and. .not.OFE_first) Return ! Vemb on disk



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
      Call mma_allocate(F_DFT,nh1,nD,Label='F_DFT')
      Call mma_allocate(D_DS,nh1,nD,Label='D_DS')
      Vxc_ref(1)=Zero
      Vxc_ref(2)=Zero
*
*---- Get the density matrix of the environment (rho_B)
*
      Call Get_iScalar('Multiplicity',kSpin)
      Call Get_D1ao(D_DS(:,1),nh1)
*     Call RecPrt('D1ao',' ',D_DS(:,1),nh1,1)
*
*---- Get the spin density matrix of the environment
*
      If (kSpin.ne.1) Then
         Call Get_D1Sao(D_DS(:,2),nh1)
*        Call RecPrt('D1Sao',' ',D_DS(:,2),nh1,1)
      End If
*
*---- Compute alpha and beta density matrices of the environment
*
      nFckDim=2
      If (kSpin.eq.1) Then
         call dscal_(nh1,Half,D_DS(:,1),1)
         call dcopy_(nh1,D_DS(:,1),1,D_DS(:,2),1)
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
*      Call RecPrt('Da',' ',D_DS(:,1),nh1,1)
*      Call RecPrt('Db',' ',D_DS(:,2),nh1,1)
      End If
*
*      If (OFE_first) Then
*---AZECH 10/2015
*   kinetic part of E_xct, Subsys B
      Func_B_TF = 0.0d0
      Call wrap_DrvNQ('TF_only',F_DFT(:,1:nFckDim),nFckDim,Func_B_TF,
     &                   D_DS(:,1:nFckDim),nh1,nFckDim,
     &                   Do_Grad,
     &                   Grad,nGrad,DFTFOCK)
*---
      If (OFE_first) Then

         Call wrap_DrvNQ(KSDFT,F_DFT(:,1:nFckDim),nFckDim,Func_B,
     &                   D_DS(:,1:nFckDim),nh1,nFckDim,
     &                   Do_Grad,
     &                   Grad,nGrad,DFTFOCK)

#ifdef _NOT_USED_
         If (KSDFT(1:4).eq.'NDSD') Then
            l_NDSD=nFckDim*nh1
            Call GetMem('NDSD','Allo','Real',ip_NDSD,l_NDSD)
            call dcopy_(l_NDSD,F_DFT(:,1:nFckDim),1,Work(ip_NDSD),1)
            KSDFT(1:4)='LDTF' !set to Thomas-Fermi for subsequent calls
         EndIf
#endif

      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*     Setup of density matrices for subsys A                           *
*                                                                      *
************************************************************************
      Call NameRun(NamRfil)    ! switch back RUNFILE name
*
      If (is_rhoA_on_file) Call NameRun('PRERFIL')
*---- Get the density matrix for rho_A
*
      Call Get_D1ao(D_DS(:,3),nh1)
*     Call RecPrt('D1ao',' ',D_DS(:,3),nh1,1)
*
      Call Get_iScalar('Multiplicity',iSpin)
      If (iSpin.eq.1 .and. kSpin.ne.1 .and. OFE_first) Then
         Call WarningMessage(0,
     &     ' Non-singlet environment perturbation on singlet state!'//
     &     '  Spin-components of the OFE potential will be averaged. ' )
      EndIf
*
*---- Get the spin density matrix of A
*
      If (iSpin.ne.1) Then
         Call Get_D1Sao(D_DS(:,4),nh1)
*        Call RecPrt('D1Sao',' ',D_DS(:,4),nh1,1)
      End If
*
*---- Compute alpha and beta density matrices of subsystem A
*
      nFckDim=2
      If (iSpin.eq.1) Then
         call dscal_(nh1,Half,D_DS(:,3),1)
         call dcopy_(nh1,D_DS(:,3),1,D_DS(:,4),1)
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
*      Call RecPrt('Da',' ',D_DS(:,3),nh1,1)
*      Call RecPrt('Db',' ',D_DS(:,4),nh1,1)
      End If
*
*---AZECH 10/2015
*   kinetic part of E_xct, Subsys A
      Call wrap_DrvNQ('TF_only',F_DFT(:,3:nFckDim+2),nFckDim,Func_A_TF,
     &                D_DS(:,3:nFckDim+2),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)
*---
      Call wrap_DrvNQ(KSDFT,F_DFT(:,3:nFckDim+2),nFckDim,Func_A,
     &                D_DS(:,3:nFckDim+2),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)
*
*  Fraction of correlation potential from A (cases: HF or Trunc. CI)
*
      If (dFMD.gt.0.0d0) Then
*
         Call mma_Allocate(Fcorr,nh1,nFckDim,Label='Fcorr')
*
         Call cwrap_DrvNQ(KSDFT,F_DFT(:,3:nFckDim+2),nFckDim,Ec_A,
     &                    D_DS(:,3:nFckDim+2),nh1,nFckDim,
     &                    Do_Grad,
     &                    Grad,nGrad,DFTFOCK,Fcorr(:,1:nFckDim))
      End If
*
*
************************************************************************
*                                                                      *
*     Calculation on the supermolecule                                 *
*                                                                      *
************************************************************************
      nFckDim=2
      If (iSpin.eq.1 .and. kSpin.eq.1) Then
         nFckDim=1
         Call daxpy_(nh1,One,D_DS(:,3),1,D_DS(:,1),1)
      Else
         Call daxpy_(nh1,One,D_DS(:,3),1,D_DS(:,1),1)
         Call daxpy_(nh1,One,D_DS(:,4),1,D_DS(:,2),1)
      EndIf
#ifdef _NOT_USED_
*---AZECH 10/2015
*   kinetic part of E_xct, Subsys A+B
*   temporarily turned off to clean output
      If (.False.) Then
       Func_AB_TF = 0.0d0
       Call wrap_DrvNQ('TF_only',F_DFT(:,1:nFckDim),nFckDim,Func_AB_TF,
     &                D_DS(:,1:nFckDim),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)
       TF_NAD = Func_AB_TF - Func_A_TF - Func_B_TF
       Write(6,*) 'kinetic part of E_xc,T (Thomas-Fermi ONLY)'
       Write(6,'(A,F19.10)') 'Ts(A+B): ', Func_AB_TF
       Write(6,'(A,F19.10)') 'Ts(A):   ', Func_A_TF
       Write(6,'(A,F19.10)') 'Ts(B):   ', Func_B_TF
       Write(6,'(A,F19.10)') '-------------------'
       Write(6,'(A,F19.10)') 'Ts_NAD:  ', TF_NAD
*   calculate v_T, Subsys A+B
       Xint_Ts_AB=dDot_(nh1,F_DFT(:,1),1,D_DS(:,3),1)
       Xint_Ts_NAD = Xint_Ts_AB - Xint_Ts_A
*     scale by 2 because wrapper only handles spin-densities
       Xint_Ts_NAD = Two*Xint_Ts_NAD
       Write(6,*) 'integrated v_Ts_NAD (Thomas-Fermi) with rhoA current'
       Write(6,'(A,F19.10)') 'Ts(A+B)_integral: ', Xint_Ts_AB
       Write(6,'(A,F19.10)') 'Ts(A)_integral:   ', Xint_Ts_A
       Write(6,'(A,F19.10)') '-------------------'
       Write(6,'(A,F19.10)') 'Ts_NAD_integral:  ', Xint_Ts_NAD
      EndIf
#endif
*---
      Call wrap_DrvNQ(KSDFT,F_DFT(:,1:nFckDim),nFckDim,Func_AB,
     &                D_DS(:,1:nFckDim),nh1,nFckDim,
     &                Do_Grad,
     &                Grad,nGrad,DFTFOCK)

      Energy_NAD = Func_AB - Func_A - Func_B
*
*---AZECH 10/2015
*   exchange-correlation part of E_xct, Subsys A+B
*   temporarily turned off to clean output
c      Write(6,*) 'E_xc_NAD (determined with Thomas-Fermi)'
c      Func_xc_NAD = Energy_NAD - TF_NAD
c      Write(6,'(A,F19.10)') 'E_xc_NAD: ', Func_xc_NAD
*---
      If (dFMD.gt.0.0d0) Then
         Call Get_electrons(xElAB)
         Fakt_ = -1.0d0*Xlambda(abs(Energy_NAD)/xElAB,Xsigma)
         Call daxpy_(nh1*nFckDim,Fakt_,Fcorr(:,1:nFckDim),1,
     &                                 F_DFT(:,3:nFckDim+2),1)
         Call mma_deallocate(Fcorr)
#ifdef _DEBUGPRINT_
         write(6,*) ' lambda(E_nad) = ',dFMD*Fakt_
#endif
      EndIf

*                                                                      *
************************************************************************
*                                                                      *
*  Non Additive (NAD) potential: F(AB)-F(A)
      Do i=1,nFckDim
         Call daxpy_(nh1,-One,F_DFT(:,2+i),1,F_DFT(:,i),1)
      End Do
#ifdef _NOT_USED_
*
*  NDSD potential for T_nad: add the (B)-dependent term
      iFickB=ip_NDSD
      Do i=1,nFckDim*Min(1,l_NDSD)
         Call daxpy_(nh1,One,Work(iFickB),1,F_DFT(:,i),1)
         If (kSpin.ne.1) iFickB=iFickB+nh1
      End Do
#endif
*
*     Add the Nuc Attr potential (from subsystem B) and then
*     put out the DFT Fock matrices from the (NAD) embedding potential
*     on the runfile (AUXRFIL). Note that the classical Coulomb
*     interaction potential from subsystem B is computed in the std
*     Fock matrix builders
*
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
*
      Call mma_allocate(TmpA,nh1,Label='TmpA')
      Call Get_dArray('Nuc Potential',TmpA,nh1)
*
      Fact = Two ! because Dmat has been scaled by half
      If (kSpin.ne.1) Fact=One
      Fact_=Fact
*
      V_emb=Fact*dDot_(nh1,F_DFT(:,1),1,D_DS(:,3),1)
      V_Nuc_AB=Fact*dDot_(nh1,TmpA,1,D_DS(:,3),1)
      If (kSpin.ne.1) Then
         V_emb=V_emb+Fact*dDot_(nh1,F_DFT(:,2),1,D_DS(:,4),1)
         V_Nuc_AB=V_Nuc_AB+Fact*dDot_(nh1,TmpA,1,D_DS(:,4),1)
      EndIf
*
*  Averaging the spin-components of F(AB) iff non-spol(A)//spol(B)
      If (iSpin.eq.1 .and. kSpin.ne.1) Then
         Do i=1,nh1
            tmp=Half*(F_DFT(i,1)+F_DFT(i,2))
            F_DFT(i,1)=tmp
         End Do
         nFckDim=1  ! reset stuff as if A+B had been spin compensated
         Fact=Two
      EndIf
*
      Do i=1,nFckDim
         Call daxpy_(nh1,1.0d0,TmpA,1,F_DFT(:,i),1)
         Vxc_ref(i)=Fact*dDot_(nh1,F_DFT(:,i),1,D_DS(:,i+2),1)
      End Do
*
      If(dFMD.gt.0.0d0) Call Put_dScalar('KSDFT energy',Ec_A)
      Call Put_dArray('Vxc_ref ',Vxc_ref,2)
*
      Call Put_dArray('dExcdRa',F_DFT(:,1:nFckDim),nh1*nFckDim)
      Call NameRun(NamRfil)   ! switch back RUNFILE name

      Call Get_dArray('Nuc Potential',TmpA,nh1)
      V_Nuc_BA= Fact_*( dDot_(nh1,TmpA,1,D_DS(:,1),1)
     &                 -dDot_(nh1,TmpA,1,D_DS(:,3),1))
      If (kSpin.ne.1) Then
         V_Nuc_BA=V_Nuc_BA+Fact_*( dDot_(nh1,TmpA,1,D_DS(:,2),1)
     &                            -dDot_(nh1,TmpA,1,D_DS(:,4),1) )
      EndIf
*
      Call mma_deallocate(TmpA)
*
#ifdef _DEBUGPRINT_
      If (nFckDim.eq.1) Then
         Do i=1,nh1
            Write(6,'(i4,f22.16)') i,F_DFT(i,1)
         End Do
      Else
         Do i=1,nh1
           Write(6,'(i4,3f22.16)') i,F_DFT(i,1),
     &                               F_DFT(i,2),
     &     (F_DFT(i,1)+F_DFT(i,2))/2.0d0
         End Do
      End If
      Write(6,'(a,f22.16)') ' NAD DFT Energy :',Energy_NAD
#endif
*
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(D_DS)
      Call Free_iSD()
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(H1)
         Call Unused_real_array(D)
         Call Unused_real(RepNuc)
         Call Unused_real(ExFac)
         Call Unused_real_array(D1I)
         Call Unused_real_array(D1A)
      End If
      End
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func,
     &                      D_DS,nh1,nD_DS,
     &                      Do_Grad,
     &                      Grad,nGrad,DFTFOCK)
      Implicit Real*8 (a-h,o-z)
      Character*(*) KSDFT
      Integer nh1, nFckDim, nD_DS
      Real*8 F_DFT(nh1,nFckDim), D_DS(nh1,nD_DS), Func
      Logical Do_Grad
      Real*8 Grad(nGrad)
      Character*4 DFTFOCK
#include "real.fh"
#include "nq_info.fh"
#include "debug.fh"
      External LSDA_emb,
     &         LSDA5_emb,
     &         BLYP_emb, BLYP_emb2,
     &         PBE_emb, PBE_emb2,
     &         Ts_only_emb, vW_hunter, nucatt_emb,
     &         Checker
      Logical  Do_MO,Do_TwoEl,F_nAsh

************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Dens_I          =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
*
      Call Get_iScalar('nSym',mIrrep)
      Call Get_iArray('nBas',mBas(0),mIrrep)
      Call Get_iArray('nFro',nFro(0),mIrrep)
      Call Get_iArray('nIsh',nIsh(0),mIrrep)
      Call qpg_dArray('nAsh',F_nAsh,nOrbA)
      If(.not.F_nAsh .or. nOrbA.eq.0) Then
         Call Izero(nAsh(0),mIrrep)
      Else
         Call Get_iArray('nAsh',nAsh(0),mIrrep)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/LSDA (Thomas-Fermi for KE)                                 *
*                                                                      *
       If (KSDFT.eq.'LDTF/LSDA ' .or.
     &     KSDFT.eq.'LDTF/LDA  ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=LDA_type
         Call DrvNQ(LSDA_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/LSDA5 (Thomas-Fermi for KE)                                *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/LSDA5' .or.
     &          KSDFT.eq.'LDTF/LDA5 ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=LDA_type
         Call DrvNQ(LSDA5_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/PBE   (Thomas-Fermi for KE)                                *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/PBE  ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=GGA_type
         Call DrvNQ(PBE_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      NDSD/PBE   (NDSD for KE)                                        *
*                                                                      *
       Else If (KSDFT.eq.'NDSD/PBE  ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=meta_GGA_type2
         Call DrvNQ(PBE_emb2,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/BLYP  (Thomas-Fermi for KE)                                *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/BLYP ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=GGA_type
         Call DrvNQ(BLYP_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      NDSD/BLYP  (NDSD for KE)                                        *
*                                                                      *
       Else If (KSDFT.eq.'NDSD/BLYP ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=meta_GGA_type2
         Call DrvNQ(BLYP_emb2,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      Kinetic only  (Thomas-Fermi)                                    *
*                                                                      *
       Else If (KSDFT.eq.'TF_only') Then
         !ExFac=Zero
         Functional_type=LDA_type
         Call DrvNQ(Ts_only_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      HUNTER  (von Weizsacker KE, no calc of potential)               *
*                                                                      *
       Else If (KSDFT.eq.'HUNTER') Then
         !ExFac=Zero
         Functional_type=GGA_type
         Call DrvNQ(vW_hunter,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      NUCATT                                                          *
*                                                                      *
       Else If (KSDFT.eq.'NUCATT_EMB') Then
         !ExFac=Zero
         Functional_type=LDA_type
         Call DrvNQ(nucatt_emb,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Checker
      Else If (KSDFT.eq.'CHECKER') Then
         !ExFac=Zero
         Functional_type=meta_GGA_type2
         Call DrvNQ(Checker,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         lKSDFT=LEN(KSDFT)
         Call WarningMessage(2,
     &               ' Wrap_DrvNQ: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
      End If
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine cWrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func,
     &                       D_DS,nh1,nD_DS,
     &                       Do_Grad,
     &                       Grad,nGrad,DFTFOCK,F_corr)
      Implicit Real*8 (a-h,o-z)
      Character*(*) KSDFT
      Integer nh1, nFckDim, nD_DS
      Real*8 F_DFT(nh1,nFckDim), D_DS(nh1,nD_DS), Func
      Real*8 F_corr(nh1,nFckDim)
      Logical Do_Grad
      Real*8 Grad(nGrad)
      Character*4 DFTFOCK
#include "real.fh"
#include "nq_info.fh"
#include "debug.fh"
      External VWN_III_emb,
     &         VWN_V_emb,
     &         cBLYP_emb,
     &         cPBE_emb,
     &         Checker
      Logical  Do_MO,Do_TwoEl,F_nAsh

************************************************************************
*                                                                      *
      Func            =Zero
      Dens_I          =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
*
      Call Get_iScalar('nSym',mIrrep)
      Call Get_iArray('nBas',mBas(0),mIrrep)
      Call Get_iArray('nFro',nFro(0),mIrrep)
      Call Get_iArray('nIsh',nIsh(0),mIrrep)
      Call qpg_dArray('nAsh',F_nAsh,nOrbA)
      If(.not.F_nAsh .or. nOrbA.eq.0) Then
         Call Izero(nAsh(0),mIrrep)
      Else
         Call Get_iArray('nAsh',nAsh(0),mIrrep)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/LSDA (Fractional) correlation potential only               *
*                                                                      *
       If (KSDFT.eq.'LDTF/LSDA ' .or.
     &     KSDFT.eq.'LDTF/LDA  ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=LDA_type
         Call DrvNQ(VWN_III_emb,F_corr,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/LSDA5 (Fractional) correlation potential only              *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/LSDA5' .or.
     &          KSDFT.eq.'LDTF/LDA5 ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=LDA_type
         Call DrvNQ(VWN_V_emb,F_corr,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/PBE   (Fractional) correlation potential only              *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/PBE  ' .or.
     &          KSDFT.eq.'NDSD/PBE  ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=GGA_type
         Call DrvNQ(cPBE_emb,F_corr,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LDTF/BLYP  (Fractional) correlation potential only              *
*                                                                      *
       Else If (KSDFT.eq.'LDTF/BLYP ' .or.
     &          KSDFT.eq.'NDSD/BLYP ') Then
         !ExFac=Get_ExFac(KSDFT(6:10))
         Functional_type=GGA_type
         Call DrvNQ(cBLYP_emb,F_corr,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Checker                                                          *
*                                                                      *
      Else If (KSDFT.eq.'CHECKER') Then
         !ExFac=Zero
         Functional_type=meta_GGA_type2
         Call DrvNQ(Checker,F_corr,nFckDim,Func,
     &              D_DS,nh1,nD_DS,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         lKSDFT=LEN(KSDFT)
         Call WarningMessage(2,
     &               ' cWrap_DrvNQ: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(F_DFT)
      End
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
      Real*8 Function Xlambda(omega,sigma)
      Implicit Real*8 (a-h,o-z)
      Real*8 omega, sigma

      If (sigma*omega.gt.42d0) Then
         Xlambda = 1.0d0
      Else
         Xlambda = 1.0d0 - exp(-sigma*omega)
      EndIf

      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Get_electrons(xnElect)
      Implicit Real*8 (a-h,o-z)
      Real*8 xnElect
#include "real.fh"
#include "nq_info.fh"

      xnElect = Dens_I

      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine VEMB_Exc_states(Vemb,nVemb,xKSDFT,Func_Bx)
      Implicit Real*8 (a-h,o-z)
      Real*8 Vemb(nVemb)
      Real*8 Func_Bx
      Character*(*) xKSDFT
      Character*16 MyNamRfil
#include <rasdim.fh>
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "ciinfo.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
#include <SysDef.fh>
      Real*8, Allocatable:: D1ao_b(:), F_DFT(:)
      Real*8, Allocatable:: xxCMO(:), xxOCCN(:), DState(:)
      Real*8 :: Dummy(1)=[0.0D0]
      Integer :: nDummy=1


      IAD12=IADR15(12)

      Call mma_allocate(xxCMO,NTOT2,Label='xxCMO')
      Call mma_allocate(xxOCCN,NTOT,Label='xxOCCN')
      Call mma_allocate(DState,NTOT1,Label='DState')
      Call mma_allocate(F_DFT,nVemb,Label='F_DFT')
      Call mma_allocate(D1ao_b,nVemb,Label='D1ao_b')

      DO KROOT=1,LROOTS
*
* Read natural orbitals
        If ( NAC.GT.0 ) then
          CALL DDAFILE(JOBIPH,2,xxCMO,NTOT2,IAD12)
          CALL DDAFILE(JOBIPH,2,xxOCCN,NTOT,IAD12)
        End If
* Get GS and excited state densities:
* Fill allocated mem with zeroes.
        DSTATE(:)=0.0D0

        Call DONE_RASSCF(xxCMO,xxOCCN,DState) ! computes D=CnC'
* Nonelectr. Vemb with GS and excited state density
        Vemb_Xstate=ddot_(nVemb,Vemb,1,DState,1)
*        Write(6,*) 'Kroot, Vemb_K ', KROOT, Vemb_Xstate
        Write(6,'(A,F19.10,3X,A,I3)') 'Nonelectr. Vemb w. rhoA_emb =',
     &        Vemb_Xstate,'root = ', KROOT
* E_xc,T[rhoA]
        Func_A=0.0d0
        F_DFT(:)=0.0D0
        Call dscal_(nVemb,0.5d0,DState,1)
        Call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_A,
     &                  DState,nVemb,1,
     &                  .false.,
     &                  Dummy,nDummy,'SCF ')
*        Write(6,*) 'Kroot, Func_A ', KROOT, Func_A
* E_xc,T[rhoA+rhoB]
        Call Get_NameRun(MyNamRfil) ! save current Runfile name
        Call NameRun('AUXRFIL')   ! switch RUNFILE name
        Call Get_D1ao(D1ao_b,nVemb)
        Call daxpy_(nVemb,0.5d0,D1ao_b,1,DState,1)
*
        Func_AB=0.0d0
        F_DFT(:)=0.0D0
        Call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_AB,
     &                  DState,nVemb,1,
     &                  .false.,
     &                  Dummy,nDummy,'SCF ')
*        Write(6,*) 'Kroot, Func_AB', KROOT, Func_AB
*        Write(6,*) 'Kroot, Func_Bx', KROOT, Func_Bx
* Calculate DFT NAD for all densities:
        DFT_NAD = Func_AB - Func_A - Func_Bx
        Write(6,'(A,F19.10,3X,A,I3)') 'DFT energy (NAD) =           ',
     &        DFT_NAD, 'root = ', KROOT
        Call NameRun(MyNamRfil) ! go back to MyNamRfil
      End Do
      Call mma_deallocate(D1ao_b)
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(DState)
      Call mma_deallocate(xxCMO)
      Call mma_deallocate(xxOCCN)

      Return
      End

