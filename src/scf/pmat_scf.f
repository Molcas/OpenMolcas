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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2003, Valera Veryazov                                  *
************************************************************************
      SubRoutine PMat_SCF(Dens,OneHam,TwoHam,nDT,NumDT,nXCf,FstItr,XCf,
     &                    nD,E_DFT,nE_DFT,Vxc,Fock)
************************************************************************
*                                                                      *
*     purpose: Compute two-electron part of the Fock matrix            *
*                                                                      *
*     input:                                                           *
*       Dens    : see SubRoutine DMat                                  *
*                                                                      *
*     output:                                                          *
*       TwoHam  : two-electron part of the Fock matrix constructed by  *
*                 contraction of proper density matrix difference with *
*                 two-electron integrals i) in conventional way ii) in *
*                 direct way                                           *
*                                                                      *
*     called from: WfCtl, Final                                        *
*                                                                      *
*     calls to: UnFold, FTwoa, Drv2El, RctFld                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: UHF - V.Veryazov, 2003                                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      External EFP_On
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "file.fh"
*
      Real*8 Dens(nDT,nD,NumDT),OneHam(nDT)
      Real*8, Dimension(:,:,:), Target:: TwoHam(nDT,nD,NumDT)
      Real*8 XCf(nXCf,nD), E_DFT(nE_DFT), Vxc(nDT,nD,NumDT)
      Real*8 Fock(nDT,nD)
      Logical FstItr, NoCoul
*
      Integer ALGO,NSCREEN
      Logical REORD,DECO
      Real*8  dmpk,dFKmat
      Common /CHOSCF / REORD,DECO,dmpk,dFKmat,ALGO,NSCREEN
*
      Logical Do_OFemb, KEonly, OFE_first, Found, EFP_On
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
*
*---- Define local variables
      Logical First, NonEq, ltmp1, ltmp2, Do_DFT
      Logical Do_ESPF
      Character NamRFil*16
      Data First /.true./
      Save First
      Real*8, Dimension(:), Allocatable:: RFfld, D
      Real*8, Dimension(:,:), Allocatable:: DnsS, Temp
      Real*8, Dimension(:,:), Allocatable, Target:: Aux
      Real*8, Dimension(:,:), Pointer:: pTwoHam
#include "SysDef.fh"

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      If (PmTime) Call CWTime(xCPM1,xWPM1)
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _DEBUG_
#ifdef _DEBUG_
      Call qEnter('PMat')
      Call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
      Call NrmClc(Vxc   (1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
      Call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
      Call NrmClc(Vxc   (1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
#endif
*
      iter_d=iter-iter0
*
* --- Copy the (abs.) value of the Max Offdiag Fmat to a Common Block
* --- Used in the LK Cholesky algorithm
      dFKmat = abs(FMOmax)
*
*---- Add contribution due to external potential
*
      Call DCopy_(nBT*nD,Zero,0,TwoHam(1,1,iPsLst),1)
      iSpin=1
      If (iUHF.eq.1) iSpin=2
      Call Put_iScalar('Multiplicity',iSpin)
*                                                                      *
************************************************************************
*                                                                      *
      Call DecideOnESPF(Do_ESPF)
      If ( Do_ESPF .or. lRF .or. KSDFT.ne.'SCF' .or. Do_OFemb .or.
     &     EFP_On()) Then
*
*------- Observe that this call always has to be prior to the calls
*        to Drv2El_dScf and/or FTwoa. This since DrvXV will redefine
*        ExFac!!!!!
*
*        Note, the linear (Oneham) and bilinear (TwoHam) contributions
*        can be computed with partial densities, as supplied with the
*        arguments to the routine. For the DFT contributions, however,
*        not being linear or bilinears, the total density is read from
*        the runfile (as put there by dmat).
*
         iCharge=Int(Tot_Charge)
         NonEq=.False.
         Do_DFT=.True.
         iDumm=1
         ltmp1=iter_d.eq.1
         ltmp2=iter_d.ne.1
         If (iUHF.eq.0) Then
            Call DrvXV(OneHam,TwoHam(1,1,iPsLst),Dens(1,1,iPsLst),
     &                  PotNuc,nBT,ltmp1,ltmp2,NonEq,
     &                  lRF,KSDFT,ExFac,iCharge,iSpin,
     &                  Dumm0,Dumm1,iDumm,'SCF ',Do_DFT)
         Else
            Call mma_allocate(D,nBT,Label='D')
            call dcopy_(nBT,Dens(1,1,iPsLst),1,D,1)
            Call DaXpY_(nBT,One,Dens(1,2,iPsLst),1,D,1)
            Call DrvXV(OneHam,TwoHam(1,1,iPsLst),D,
     &                 PotNuc,nBT,ltmp1,ltmp2,NonEq,
     &                 lRF,KSDFT,ExFac,iCharge,iSpin,
     &                 Dumm0,Dumm1,iDumm,'SCF ',Do_DFT)
            Call mma_deallocate(D)
            call dcopy_(nBT,TwoHam(1,1,iPsLst),1,TwoHam(1,2,iPsLst),1)
            If (MxConstr.gt.0 .and. klockan.eq.1) Then
               Call SetUp_iSD()
               Call Get_Enondyn_dft(nBT,Dumm1,iDumm,'SCF ')
               Call Free_iSD()
               klockan=24
            EndIf
         End If
*
*        Pick up the integrated energy contribution of the external
*        potential to the total energy.
*
         Call Peek_dScalar('KSDFT energy',E_DFT(iter_d))
*
*        Pick up the contribution to the Fock matrix due to the
*        external field. Note that for some external field the
*        potential is neither linear nor bi-linear.
*
         If (KSDFT.ne.'SCF') Then
            Call Get_dExcdRa(ipVxc,nVxc)
            Call DCopy_(nVxc,Work(ipVxc),1,Vxc(1,1,iPsLst),1)
            Call Free_Work(ipVxc)
         Else
            Call FZero(Vxc(1,1,iPsLst),nBT*nD)
         End If
*
         If (Do_OFemb) Then
            Call Get_NameRun(NamRfil) ! save the old RUNFILE name
            Call NameRun('AUXRFIL')   ! switch the RUNFILE name
            Call Get_dExcdRa(ipVemb,nVemb)
            Call DaXpY_(nDT*nD,One,Work(ipVemb),1,Vxc(1,1,iPsLst),1)
            Call Free_Work(ipVemb)
            Call NameRun(NamRfil)   ! switch back RUNFILE name
         End If
#ifdef _DEBUG_
         Call NrmClc(Vxc   (1,1,iPsLst),nDT*nD,'PMat','Optimal V ')
#endif
*
      Else If ( RFpert.and.First ) Then
*
         If (iUHF.eq.1) Then
            write(6,*) ' UHF+RF: Not implemented'
            call Abend()
         End If
         Call mma_allocate(RFfld,nBT,Label='RFfld')
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call Get_dScalar('RF Self Energy',ERFself)
         Call Get_dArray('Reaction field',RFfld,nBT)
         If (Found) Call NameRun('RUNFILE')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nBT,1.0d0,RFfld,1,OneHam,1)
         Do iD = 1, nD
            Call DCopy_(nBT,OneHam,1,Fock(1,iD),1)
         End Do
         Call mma_deallocate(RFfld)
*
      Else
*
         Call FZero(Vxc(1,1,iPsLst),nBT*nD)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      First=.false.
*
*---- Compute the two-electron contribution to the Fock matrix
*
      Backup_ExFac=ExFac
      If (NoExchange) Then
         ExFac=0.0d0
      End If
      nT = 1
      If (nD.eq.2) nT=3
      Call mma_allocate(Temp,nBT,nT,Label='Temp')
      Call FZero(Temp,nBT*nT)
      If (PmTime) Call CWTime(tCF2,tWF2)
      If (DSCF) Then
         If (iUHF.eq.0) Then
            NoCoul=.False.
            Call Drv2El_dscf(Dens(1,1,iPsLst),Temp(1,1),nBT,
     &                       Max(nDisc*1024,nCore),Thize,PreSch,FstItr,
     &                       NoCoul,ExFac)
         Else
*
*           Compute the Coulomb potential for the total density and
*           exchange of alpha and beta, respectively. Add together
*           to get the correct contributions to the alpha and beta
*           Fock matrices.
*
*           Set exchange factor to zero and compute only Coulomb
*           for the total electron density.
*
            NoCoul=.False.
            Call DCopy_(nBT,Dens(1,1,iPsLst),1,Temp(1,2),1)
            Call DaXpY_(nBT,1.0D0,Dens(1,2,iPsLst),1,Temp(1,2),1)
*
            Call Drv2El_dscf(Temp(1,2),Temp(1,3),nBT,
     &                       Max(nDisc*1024,nCore),Thize,PreSch,FstItr,
     &                       NoCoul,0.0D0)
*
*           alpha exchange
            NoCoul=.TRUE.
            Call FZero(Temp(1,2),nBT)
            Call Drv2El_dscf(Dens(1,1,iPsLst),Temp(1,1),nBT,
     &                       Max(nDisc*1024,nCore),Thize,PreSch,FstItr,
     &                       NoCoul,ExFac)
            Call DScal_(nBT,2.0D0,Temp(1,1),1)
*
*           beta exchange
            Call Drv2El_dscf(Dens(1,2,iPsLst),Temp(1,2),nBT,
     &                       Max(nDisc*1024,nCore),Thize,PreSch,FstItr,
     &                       NoCoul,ExFac)
            Call DScal_(nBT,2.0D0,Temp(1,2),1)
*
*           Add together J and K contributions to form the correct
*           alpha and beta Fock matrices.
*
            Call DaXpY_(nBT,1.0D0,Temp(1,3),1,Temp(1,1),1)
            Call DaXpY_(nBT,1.0D0,Temp(1,3),1,Temp(1,2),1)
         End If
      Else
*------- Allocate memory for squared density matrix
         Call mma_allocate(DnsS,nBB,nD,Label='DnsS')
*------- Expand the 1-density matrix
         Do iD = 1, nD
            Call Unfold(Dens(1,iD,iPsLst),nBT,DnsS(1,iD),nBB,nSym,nBas)
         End Do
         If (iUHF.eq.0) Then
            Call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,
     &                     Dens(1,1,iPsLst),DnsS(1,1),Temp(1,1),
     &                     nBT,ExFac,nBB,MaxBas,iUHF,
     &                     Dummy,
     &                     Dummy,Dummy,nOcc(1,1),idummy,
     &                     iDummy_run)
         Else
            Call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,
     &                     Dens(1,1,iPsLst),DnsS(1,1),Temp(1,1),
     &                     nBT,ExFac,nBB,MaxBas,iUHF,
     &                     Dens(1,2,iPsLst),
     &                     DnsS(1,2),Temp(1,2),nOcc(1,1),
     &                     nOcc(1,2),iDummy_run)
         End If
*
*------- Deallocate memory for squared density matrix
         Call mma_deallocate(DnsS)
      End If
      If (PmTime) Then
         Call CWTime(tCF2_1,tWF2_1)
         tCF2=tCF2_1-tCF2
         tWF2=tWF2_1-tWF2
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Add on FMM contributions to Fock matrix, only works for HF
*
      If (DoFMM) Call FMMFck(Dens(1,1,iPsLst),Temp(1,1),nDens)
*                                                                      *
************************************************************************
*                                                                      *
      Call DaXpY_(nBT*nD,One,Temp,1,TwoHam(1,1,iPsLst),1)
#ifdef _DEBUG_
      Call NrmClc(Temp,nBT*nD,'PMat_SCF','Temp')
      Call NrmClc(TwoHam(1,1,iPsLst),nBT*nD,'PMat_SCF','T in iPsLst')
#endif
      Call mma_deallocate(Temp)
*                                                                      *
************************************************************************
*                                                                      *
*     Now compute the total two-electron contribution
*
*     Generate the two-electron contibution corresponding to the total
*     density.
*
      If(MiniDn.and. Max(0,nIter(nIterP)-1).gt.0) Then
*
*        Minimized density option
*
*        D(k+1) = Sum_i_k C_i D_i + delta(k+1)
*
*        G(D(k+1)) = G(delta(k+1)) + Sum_i_k C_i G(D_i)
*
         Call mma_allocate(Aux,nDT,nD,Label='Aux')
         Do iMat = 1, iter_d-1
*
            tmp = 0.0D0
            Do iD = 1, nD
               tmp = tmp + Abs(XCf(iMat,iD))
            End Do
            If (tmp.eq.0.0D0) Cycle

            iM = MapDns(iMat)
            If (iM.lt.0) Then
               Call RWDTG(-iM,Aux,nDT*nD,'R','TWOHAM',iDisk,MxDDsk)
               pTwoHam => Aux
            Else
               pTwoHam => TwoHam(1:nDT,1:nD,iM)
            End If
*
            Do iD = 1, nD
               If (Xcf(iMat,iD).eq.0.0D0) Cycle
               Call DaXpY_(nBT,Xcf(iMat,iD),pTwoHam(1,iD),1,
     &                              TwoHam(1,iD,iPsLst),1)
            End Do
*
            Nullify(pTwoHam)
*
         End Do
         Call mma_deallocate(Aux)
*
      Else If (.not.DDnOFF) Then
*
*        Normal density difference
*
*        G(D(k+1)) = G(D(k)) + G(D(k+1)-D(k))
*
         Call DaXpY_(nBT*nD,One,TwoHam(1,1,nDens ),1,
     &                          TwoHam(1,1,iPsLst),1)
*
      End If
*
*     Restore the total density in position iPsLst
*
      Call DCopy_(nBT*nD,Dens  (1,1,nDens ),1,Dens  (1,1,iPsLst),1)
      Call DCopy_(nBT*nD,TwoHam(1,1,iPsLst),1,TwoHam(1,1,nDens ),1)
      Call DCopy_(nBT*nD,Vxc   (1,1,iPsLst),1,Vxc   (1,1,nDens ),1)
*
      ! Restore ExFac (if it was changed)
      If (NoExchange) Then
         ExFac=Backup_ExFac
      End If
*                                                                      *
************************************************************************
*                                                                      *
      TNorm = DDot_(nBT*nD,TwoHam(1,1,iPsLst),1,TwoHam(1,1,iPsLst),1)
     &      / DBLE(nD)
*
*
*define _DEBUG_
#ifdef _DEBUG_
      Call NrmClc(Dens  (1,1,iPsLst),nBT*nD,'PMat  ','D iPsLst  ')
      Call NrmClc(Dens  (1,1,nDens), nBT*nD,'PMat  ','D nDens   ')
      Call NrmClc(TwoHam(1,1,iPsLst),nBT*nD,'PMat  ','T iPsLst  ')
      Call NrmClc(TwoHam(1,1,nDens), nBT*nD,'PMat  ','T nDens   ')
      Call NrmClc(Vxc   (1,1,iPsLst),nBT*nD,'PMat  ','V iPsLst  ')
      Call NrmClc(Vxc   (1,1,nDens), nBT*nD,'PMat  ','V nDens   ')
*
      Call qExit('PMat')
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 5) = TimFld( 5) + (Cpu2 - Cpu1)
      If (PmTime) Then
         Call CWTime(xCPM2,xWPM2)
         xCPM = xCPM2 - xCPM1
         xWPM = xWPM2 - xWPM1
         Write(6,'(1X,A,F15.2,A,A,F15.2,A,/,1X,A,F15.2,A,A,F15.2,A)')
     &   '>>> PMat_SCF: CPU  time:',xCPM,' seconds ',
     &   ' (2-el contributions: ',tCF2,' seconds) <<<',
     &   '>>> PMat_SCF: Wall time:',xWPM,' seconds ',
     &   ' (2-el contributions: ',tWF2,' seconds) <<<'
         Call xFlush(6)
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
