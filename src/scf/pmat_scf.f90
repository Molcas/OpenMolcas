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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!***********************************************************************
      SubRoutine PMat_SCF(FstItr,XCf,nXCF,nD)
!***********************************************************************
!                                                                      *
!     purpose: Compute two-electron part of the Fock matrix            *
!                                                                      *
!     output:                                                          *
!       TwoHam  : two-electron part of the Fock matrix constructed by  *
!                 contraction of proper density matrix difference with *
!                 two-electron integrals i) in conventional way ii) in *
!                 direct way                                           *
!                                                                      *
!***********************************************************************
      use OFembed, only: Do_OFemb
      use InfSCF, Only: pMTime, nBT, ipsLst, KSDFT, Iter, MxConstr, Klockan, RFPert, PotNuc, exFac,          &
                        NoExchange, DSCF, nDIsc, nOcc, DoFMM, MiniDn, nDens, tNorm, DDnOff, FMOMax, iDisk,         &
                        iDummy_Run, MapDns, MaxBas, nBas, nBB, nCore, nIter, nIterP, nSkip, nSym, PreSch, Thize,   &
                        TimFld, Tot_Charge
      use ChoSCF, only: dfkMat, Algo
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One, Two
      use RICD_Info, only: Do_DCCD
      use SCF_Arrays, only: Dens, OneHam, TwoHam, Vxc, Fock=>FockAO, EDFT
      use Int_Options, only: FckNoClmb, Exfac_Int=>ExFac, Thize_Int=>Thize, PreSch_Int=>PreSch
      use rctfld_module, only: lRF
      Implicit None
      External EFP_On
!
      Integer nD, nXCf
      Real*8 XCf(nXCf,nD)
!
!---- Define local variables
      Integer nDT, iSpin, iCharge, iDumm, nVxc, iD, nT, iMat, iM
      Integer Algo_Save
      Logical FstItr, Found, EFP_On
      Logical, Save :: First=.True.
      Logical NonEq, ltmp1, ltmp2, Do_DFT, Do_ESPF
      Real*8 ERFSelf, Backup_ExFac, TCF2, TCF2_1, TWF2, TWF2_1
      Real*8 Tmp, CPU1, CPU2, XCPM, XCPM1, XCPM2
      Real*8                  XWPM, XwPM1, XwPM2
      Real*8 Tim1, Tim2, Tim3
      Real*8, Dimension(:), Allocatable:: RFfld, D
      Real*8, Dimension(:,:), Allocatable:: DnsS
      Real*8, Allocatable, Target:: Temp(:,:)
      Real*8, Dimension(:,:), Allocatable, Target:: Aux
      Real*8, Dimension(:,:), Pointer:: pTwoHam
      Real*8, Allocatable :: tVxc(:)
      Real*8, External :: DDot_
      Real*8 Dummy(1),Dumm0(1),Dumm1(1)
      Real*8, Allocatable:: Save(:,:)
#include "SysDef.fh"
!
      Interface
        SubRoutine Drv2El_dscf(Dens,TwoHam,nDens,nDisc,FstItr)
        Integer nDens, nDisc
        Real*8, Target:: Dens(nDens), TwoHam(nDens)
        Logical FstItr
        End Subroutine Drv2El_dscf
      End Interface

      nDT = Size(OneHam)
      If (PmTime) Call CWTime(xCPM1,xWPM1)
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
#ifdef _DEBUGPRINT_
      Call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
      Call NrmClc(Vxc   (1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
      Call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
      Call NrmClc(Vxc   (1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
#endif
!
! --- Copy the (abs.) value of the Max Offdiag Fmat to a Common Block
! --- Used in the LK Cholesky algorithm
      dFKmat = abs(FMOmax)
!
!---- Add contribution due to external potential
!
      Call DCopy_(nBT*nD,[Zero],0,TwoHam(1,1,iPsLst),1)
      iSpin=1
      If (nD==2) iSpin=2
      Call Put_iScalar('Multiplicity',iSpin)
!                                                                      *
!***********************************************************************
!                                                                      *
      Call DecideOnESPF(Do_ESPF)
      If ( Do_ESPF .or. lRF .or. KSDFT.ne.'SCF' .or. Do_OFemb .or. EFP_On()) Then
!
!------- Observe that this call always has to be prior to the calls
!        to Drv2El_dScf and/or FTwoa. This since DrvXV will redefine
!        ExFac!!!!!
!
!        Note, the linear (Oneham) and bilinear (TwoHam) contributions
!        can be computed with partial densities, as supplied with the
!        arguments to the routine. For the DFT contributions, however,
!        not being linear or bilinears, the total density is read from
!        the runfile (as put there by dmat).
!
         iCharge=Int(Tot_Charge)
         NonEq=.False.
         Do_DFT=.True.
         iDumm=1
         ltmp1=iter.eq.1
         ltmp2=iter.ne.1
         If (nD==1) Then
            Call DrvXV(OneHam,TwoHam(1,1,iPsLst),Dens(1,1,iPsLst),PotNuc,nBT,ltmp1,ltmp2,NonEq,    &
                        lRF,KSDFT,ExFac,iCharge,iSpin,Dumm0,Dumm1,iDumm,'SCF ',Do_DFT)
         Else
            Call mma_allocate(D,nBT,Label='D')
            call dcopy_(nBT,Dens(1,1,iPsLst),1,D,1)
            Call DaXpY_(nBT,One,Dens(1,2,iPsLst),1,D,1)
            Call DrvXV(OneHam,TwoHam(1,1,iPsLst),D,PotNuc,nBT,ltmp1,ltmp2,NonEq,       &
                       lRF,KSDFT,ExFac,iCharge,iSpin,Dumm0,Dumm1,iDumm,'SCF ',Do_DFT)
            Call mma_deallocate(D)
            call dcopy_(nBT,TwoHam(1,1,iPsLst),1,TwoHam(1,2,iPsLst),1)
            If (MxConstr.gt.0 .and. klockan.eq.1) Then
               Call SetUp_iSD()
               Call Get_Enondyn_dft(nBT,Dummy,iDumm,'SCF ')
               Call Free_iSD()
               klockan=24
            EndIf
         End If
!
!        Pick up the integrated energy contribution of the external
!        potential to the total energy.
!
         Call Peek_dScalar('KSDFT energy',EDFT(iter))
!
!        Pick up the contribution to the Fock matrix due to the
!        external field. Note that for some external field the
!        potential is neither linear nor bi-linear.
!
         If (KSDFT.ne.'SCF') Then
            nVxc=Size(Vxc,1)*Size(Vxc,2)
            Call mma_allocate(tVxc,nVxc,Label='tVxc')
            Call Get_dArray_chk('dExcdRa',tVxc,nVxc)
            Call DCopy_(nVxc,tVxc,1,Vxc(1,1,iPsLst),1)
            Call mma_deallocate(tVxc)
         Else
            Call FZero(Vxc(1,1,iPsLst),nBT*nD)
         End If
!
         If (Do_OFemb) Then
            Call NameRun('AUXRFIL') ! switch the RUNFILE name
            nVxc=Size(Vxc,1)*Size(Vxc,2)
            Call mma_allocate(tVxc,nVxc,Label='tVxc')
            Call Get_dArray_chk('dExcdRa',tVxc,nVxc)
            Call DaXpY_(nDT*nD,One,tVxc,1,Vxc(1,1,iPsLst),1)
            Call mma_deallocate(tVxc)
            Call NameRun('#Pop')    ! switch back RUNFILE name
         End If
#ifdef _DEBUGPRINT_
         Call NrmClc(Vxc   (1,1,iPsLst),nDT*nD,'PMat','Optimal V ')
#endif
!
      Else If ( RFpert.and.First ) Then
!
         If (nD==2) Then
            write(6,*) ' UHF+RF: Not implemented'
            call Abend()
         End If
         Call mma_allocate(RFfld,nBT,Label='RFfld')
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call Get_dScalar('RF Self Energy',ERFself)
         Call Get_dArray('Reaction field',RFfld,nBT)
         If (Found) Call NameRun('#Pop')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nBT,One,RFfld,1,OneHam,1)
         Do iD = 1, nD
            Call DCopy_(nBT,OneHam,1,Fock(1,iD),1)
         End Do
         Call mma_deallocate(RFfld)
!
      Else
!
         Call FZero(Vxc(1,1,iPsLst),nBT*nD)
!
      End If
!                                                                      *
!***********************************************************************

      First=.false.
!
!---- Compute the two-electron contribution to the Fock matrix
!
      Backup_ExFac=ExFac
      If (NoExchange) Then
         ExFac=Zero
      End If

      nT = 1 + (nD-1)*2

      Call mma_allocate(Temp,nBT,nT,Label='Temp')
      Temp(:,:)=Zero

      If (PmTime) Call CWTime(tCF2,tWF2)

      If (DSCF.and..Not.Do_DCCD) Then

         Call Drv2El_dscf_Front_End(Temp,Size(Temp,1),Size(Temp,2))

      Else   ! RICD/Cholesky option

!------- Allocate memory for squared density matrix
         Call mma_allocate(DnsS,nBB,nD,Label='DnsS')

!------- Expand the 1-density matrix
         Do iD = 1, nD
            Call Unfold(Dens(1,iD,iPsLst),nBT,DnsS(1,iD),nBB,nSym,nBas)
         End Do

         Call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,Dens(:,:,iPsLst),DnsS(:,:),Temp(:,:),        &
                              nBT,ExFac,nBB,MaxBas,nD,nOcc(:,:),Size(nOcc,1),iDummy_run)

         If (Do_DCCD) Then
            Call mma_Allocate(Save,Size(Temp,1),Size(Temp,2),Label='Save')

            Save(:,:)=Zero
            Call Drv2El_dscf_Front_End(Save,Size(Temp,1),Size(Temp,2))
            Temp(:,:) = Temp(:,:) + Save(:,:)

            Algo_save=Algo
            Algo=0
            Save(:,:)=Zero
            Call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,Dens(:,:,iPsLst),DnsS(:,:),Save(:,:),      &
                                 nBT,ExFac,nBB,MaxBas,nD,nOcc(:,:),Size(nOcc,1),iDummy_run)
            Temp(:,:) = Temp(:,:) - Save(:,:)
            Algo=Algo_save
            Call mma_deAllocate(Save)
         End If
!
!------- Deallocate memory for squared density matrix
         Call mma_deallocate(DnsS)

      End If

      If (PmTime) Then
         Call CWTime(tCF2_1,tWF2_1)
         tCF2=tCF2_1-tCF2
         tWF2=tWF2_1-tWF2
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Add on FMM contributions to Fock matrix, only works for HF
!
      If (DoFMM) Call FMMFck(Dens(1,1,iPsLst),Temp(1,1),nDens)
!                                                                      *
!***********************************************************************
!                                                                      *
      Call DaXpY_(nBT*nD,One,Temp,1,TwoHam(1,1,iPsLst),1)
#ifdef _DEBUGPRINT_
      Call NrmClc(Temp,nBT*nD,'PMat_SCF','Temp')
      Call NrmClc(TwoHam(1,1,iPsLst),nBT*nD,'PMat_SCF','T in iPsLst')
#endif
      Call mma_deallocate(Temp)

!                                                                      *
!***********************************************************************
!                                                                      *
!     Now compute the total two-electron contribution
!
!     Generate the two-electron contribution corresponding to the total
!     density.
!
      If(MiniDn.and. Max(0,nIter(nIterP)-1).gt.0) Then
!
!        Minimized density option
!
!        D(k+1) = Sum_i_k C_i D_i + delta(k+1)
!
!        G(D(k+1)) = G(delta(k+1)) + Sum_i_k C_i G(D_i)
!
         Call mma_allocate(Aux,nDT,nD,Label='Aux')
         Do iMat = 1, iter-1
!
            tmp = Zero
            Do iD = 1, nD
               tmp = tmp + Abs(XCf(iMat,iD))
            End Do
            If (tmp.eq.Zero) Cycle

            iM = MapDns(iMat)
            If (iM.lt.0) Then
               Call RWDTG(-iM,Aux,nDT*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
               pTwoHam => Aux
            Else
               pTwoHam => TwoHam(1:nDT,1:nD,iM)
            End If
!
            Do iD = 1, nD
               If (Xcf(iMat,iD).eq.Zero) Cycle
               Call DaXpY_(nBT,Xcf(iMat,iD),pTwoHam(:,iD),1,TwoHam(1,iD,iPsLst),1)
            End Do
!
            Nullify(pTwoHam)
!
         End Do
         Call mma_deallocate(Aux)
!
      Else If (.not.DDnOFF) Then
!
!        Normal density difference
!
!        G(D(k+1)) = G(D(k)) + G(D(k+1)-D(k))
!
         Call DaXpY_(nBT*nD,One,TwoHam(1,1,nDens ),1,TwoHam(1,1,iPsLst),1)
!
      End If
!
!     Restore the total density in position iPsLst
!
      Call DCopy_(nBT*nD,Dens  (1,1,nDens ),1,Dens  (1,1,iPsLst),1)
      Call DCopy_(nBT*nD,TwoHam(1,1,iPsLst),1,TwoHam(1,1,nDens ),1)
      Call DCopy_(nBT*nD,Vxc   (1,1,iPsLst),1,Vxc   (1,1,nDens ),1)
!
      ! Restore ExFac (if it was changed)
      If (NoExchange) Then
         ExFac=Backup_ExFac
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
      TNorm = DDot_(nBT*nD,TwoHam(1,1,iPsLst),1,TwoHam(1,1,iPsLst),1) / DBLE(nD)
!
!
#ifdef _DEBUGPRINT_
      Call NrmClc(Dens  (1,1,iPsLst),nBT*nD,'PMat  ','D iPsLst  ')
      Call NrmClc(Dens  (1,1,nDens), nBT*nD,'PMat  ','D nDens   ')
      Call NrmClc(TwoHam(1,1,iPsLst),nBT*nD,'PMat  ','T iPsLst  ')
      Call NrmClc(TwoHam(1,1,nDens), nBT*nD,'PMat  ','T nDens   ')
      Call NrmClc(Vxc   (1,1,iPsLst),nBT*nD,'PMat  ','V iPsLst  ')
      Call NrmClc(Vxc   (1,1,nDens), nBT*nD,'PMat  ','V nDens   ')
!
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 5) = TimFld( 5) + (Cpu2 - Cpu1)
      If (PmTime) Then
         Call CWTime(xCPM2,xWPM2)
         xCPM = xCPM2 - xCPM1
         xWPM = xWPM2 - xWPM1
         Write(6,'(1X,A,F15.2,A,A,F15.2,A,/,1X,A,F15.2,A,A,F15.2,A)')         &
         '>>> PMat_SCF: CPU  time:',xCPM,' seconds ',                         &
         ' (2-el contributions: ',tCF2,' seconds) <<<',                       &
         '>>> PMat_SCF: Wall time:',xWPM,' seconds ',                         &
         ' (2-el contributions: ',tWF2,' seconds) <<<'
         Call xFlush(6)
      End If

      Contains
      Subroutine Drv2El_dscf_Front_End(Temp,n1,n2)
      Integer n1, n2
      Real*8 Temp(n1,n2)

! while the Drv2El_dscf can't handle UHF in a trivial way this interface has
! to be used.

      ExFac_Int=ExFac
      Thize_Int=Thize
      PreSch_Int=PreSch
      If (n2==1) Then
         FckNoClmb=.False.
         ExFac_Int=ExFac
         Call Drv2El_dscf(Dens(1,1,iPsLst),Temp(1,1),nBT,0,FstItr)
      Else
!
!        Compute the Coulomb potential for the total density and
!        exchange of alpha and beta, respectively. Add together
!        to get the correct contributions to the alpha and beta
!        Fock matrices.
!
!        Set exchange factor to zero and compute only Coulomb
!        for the total electron density.
!
         FckNoClmb=.False.
         Temp(:,2)=Dens(:,1,iPsLst)+Dens(:,2,iPsLst)
!
         ExFac_Int=Zero
         Call Drv2El_dscf(Temp(1,2),Temp(1,3),nBT,Max(nDisc*1024,nCore),FstItr)
!
!        alpha exchange
         FckNoClmb=.TRUE.
         Temp(:,2)=Zero
         ExFac_Int=ExFac
         Call Drv2El_dscf(Dens(1,1,iPsLst),Temp(1,1),nBT,Max(nDisc*1024,nCore),FstItr)
         Temp(:,1)=Two*Temp(:,1)
!
!        beta exchange
         ExFac_Int=ExFac
         Call Drv2El_dscf(Dens(1,2,iPsLst),Temp(1,2),nBT,Max(nDisc*1024,nCore),FstItr)
         Temp(:,2)=Two*Temp(:,2)
!
!        Add together J and K contributions to form the correct
!        alpha and beta Fock matrices.
!
         Temp(:,1)=Temp(:,1)+Temp(:,3)
         Temp(:,2)=Temp(:,2)+Temp(:,3)
       End If

       End Subroutine Drv2El_dscf_Front_End

      End subroutine PMat_SCF
