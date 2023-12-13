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
! Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
!               1990, IBM                                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Drv2El_dscf(Dens,TwoHam,nDens,nDisc,FstItr)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will control the type of the two-electron integral, e.g.    *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
!          Dens is the folded lower triangular of the 1st order        *
!               density matrix.                                        *
!          Twoham is the lower triangular of the two-electron contri-  *
!               bution to the Fock matrix.                             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for direct SCF. January '93                     *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified by M.Schuetz @teokem.lu.se :                    *
!             parallel region split off in drvtwo.f, April '95         *
!             Modified by R. Lindh  @teokem.lu.se :                    *
!             total repacking of code September '96                    *
!***********************************************************************
      use IOBUF, only: lBuf
      use Gateway_Info, only: ThrInt, CutInt
      use RICD_Info, only: Do_DCCD
      use iSD_data, only: iSD
      use Integral_Interfaces, only: DeDe_SCF, No_routine,
     &                               Int_PostProcess
      use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
      use Int_Options, only: Exfac, Thize, W2Disc
      use Int_Options, only: Disc_Mx, Disc, Count=>Quad_ijkl
      use Constants, only: Zero, One, Two, Three, Four, Eight
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nDens, nDisc
      Real*8, Target:: Dens(nDens), TwoHam(nDens)
      Logical FstItr
!
      External Rsv_GTList
      Integer, Parameter :: nTInt=1
      Real*8 TInt(nTInt)
      Logical Semi_Direct,Rsv_GTList, Indexation,
     &        DoGrad, Triangular
      Character(LEN=72) SLine
      Real*8, Allocatable:: TMax(:,:), DMax(:,:)
      Integer, Allocatable:: ip_ij(:,:)
      Integer iS, jS, ijS, klS, nSkal, iOpt, nIJ, kS, lS, mDens
      Real*8 ThrAO, TskHi, TskLw, P_Eff, PP_Eff, PP_Eff_Delta,
     &       PP_Count, TMax_All, S_Eff, T_Eff, ST_Eff, AInt, Dtst,
     &       TCPU1, TCPU2, TWALL1, TWALL2
!                                                                      *
!***********************************************************************
!                                                                      *
      SLine='Computing 2-electron integrals'
      Call StatusLine(' SCF:',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      Int_PostProcess => No_Routine
!                                                                      *
!***********************************************************************
!                                                                      *
!     Set variables in module Int_Options
      DoIntegrals=.False.
      DoFock=.True.
      FckNoExch=ExFac.eq.Zero
      W2Disc=.False.     ! Default value
!     Disc_Mx = file size in Real*8 128=1024/8
      Disc_Mx= DBLE(nDisc)*128.D00
!     Subtract for the last buffer
      Disc_Mx= Disc_Mx - lBuf
      Disc = Zero        ! Default value
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Set up for partial SO/AO integral storage.
!
!---- nDisc = file size in kbyte from input
      Semi_Direct = nDisc.ne.0
      If (Semi_Direct) Then
         Call Init_SemiDSCF(FstItr,Thize,Cutint)
      Endif
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Desymmetrize differential densities.
!     Observe that the desymmetrized 1st order density matrices are
!     canonical, i.e. the relative order of the indices are canonically
!     ordered.
!
      Call DeDe_SCF(Dens,TwoHam,nDens,mDens)
!                                                                      *
!***********************************************************************
!                                                                      *
      Indexation=.False.
      ThrAO=Zero           ! Do not modify CutInt
      DoGrad=.False.
!
      Call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
      TskHi=Zero
      TskLw=Zero
!                                                                      *
!***********************************************************************
!                                                                      *
!---  Compute entities for prescreening at shell level
!
      Call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
      Call Shell_MxSchwz(nSkal,TMax)
      TMax_all=Zero
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (Do_DCCD.and.iSD(10,iS)/=iSD(10,jS)) Cycle
            TMax_all=Max(TMax_all,TMax(iS,jS))
         End Do
      End Do
      Call mma_allocate(DMax,nSkal,nSkal,Label='DMax')
      Call Shell_MxDens(Dens,DMax,nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Create list of non-vanishing pairs
!
      Call mma_allocate(ip_ij,2,nSkal*(nSkal+1),Label='ip_ij')
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
         If (Do_DCCD.and.iSD(10,iS)/=iSD(10,jS)) Cycle
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nij = nij + 1
               ip_ij(1,nij)=iS
               ip_ij(2,nij)=jS
            End If
         End Do
      End Do
      P_Eff=Dble(nij)
!
      PP_Eff=P_Eff**2
      PP_Eff_delta=0.10D0*PP_Eff
      PP_Count=Zero

!                                                                      *
!***********************************************************************
!                                                                      *
!.... For distributed parallel SCF initiate (sequential code is special
!     case when the number of nodes in the mpp is 1).
!
!     1: Task list (tlist)
!     2: Private priority list (pplist)
!     3: Global task list (gtlist)
!
      If (FstItr) Then
         Triangular=.True.
         Call Init_TList(Triangular,P_Eff)
         Call Init_PPList
         Call Init_GTList
      Else
         Call ReInit_PPList(Semi_Direct)
         Call ReInit_GTList
      End If
      iOpt=0
      If (.Not.FstItr.and.Semi_direct) Then
         iOpt=2
      Endif
!
      Call CWTime(TCpu1,TWall1)
!
!     big loop over individual tasks, distributed over individual nodes

   10 Continue
!     make reservation of a task on global task list and get task range
!     in return. Function will be false if no more tasks to execute.

      If (.Not.Rsv_GTList(TskLw,TskHi,iOpt,W2Disc)) Then
         Go To 11
      Endif


      Call Mode_SemiDSCF(W2Disc)
!     Write (6,*) 'TskLw,TskHi,W2Disc=',TskLw,TskHi,W2Disc
!
!     Now do a quadruple loop over shells
!
      ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
      iS = ip_ij(1,ijS)
      jS = ip_ij(2,ijS)
      klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      kS = ip_ij(1,klS)
      lS = ip_ij(2,klS)
      Count=TskLw


      If (Count-TskHi.gt.1.0D-10) Go To 12 ! Cut off check
! What are these variables
  13  Continue
!
      S_Eff=DBLE(ijS)
      T_Eff=DBLE(klS)
      ST_Eff=S_Eff*(S_Eff-One)/2D0 + T_Eff


      If (ST_Eff.ge.PP_Count) Then
         Write (SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',
     &        ST_Eff/PP_Eff,'% done so far.'
         Call StatusLine(' Seward:',SLine)
         PP_Count = PP_Count + PP_Eff_delta
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
         Aint=TMax(iS,jS)*TMax(kS,lS)
         If (Semi_Direct) Then

!
!           No density screening in semi-direct case!
!           Cutint: Threshold for Integrals. In semi-direct case, this
!                   must be the final threshold used in the last scf
!                   iteration
!           Thrint: Threshold for Density screening. This the actual
!                   threshold
!                   for the current iteration
!
           If (AInt.lt.CutInt) Go To 14
         Else

           If(FckNoClmb) then
              Dtst=Max(DMax(is,ls)/Four,DMax(is,ks)/Four,
     &                 DMax(js,ls)/Four,DMax(js,ks)/Four)
           Else If(FckNoExch) then
              Dtst=Max(DMax(is,js),DMax(ks,ls))
           Else
              Dtst=Max(DMax(is,ls)/Four,DMax(is,ks)/Four,
     &                 DMax(js,ls)/Four,DMax(js,ks)/Four,
     &                 DMax(is,js),DMax(ks,ls))
           End If

           If (Aint*Dtst.lt.ThrInt) Then
              goto 14
           Endif

         End if
         If (Do_DCCD.and.iSD(10,iS)/=iSD(10,kS)) Go To 14
!                                                                      *
!***********************************************************************
!                                                                      *
         Call Eval_IJKL(iS,jS,kS,lS,TInt,nTInt)

 14      Continue
         Count=Count+One
         If (Count-TskHi.gt.1.0D-10) Go To 12
         klS = klS + 1
         If (klS.gt.ijS) Then
            ijS = ijS + 1
            klS = 1
         End If
         iS = ip_ij(1,ijS)
         jS = ip_ij(2,ijS)
         kS = ip_ij(1,klS)
         lS = ip_ij(2,klS)
         Go To 13
!
!     Task endpoint
!
 12   Continue
!
      If (Semi_Direct) Then
         If (W2Disc) Then
            Call Put_QLast
         Else
            Call Pos_QLast(Disc)
         End If
      End If
!
      Go To 10
 11   Continue
!     End of big task loop
      Call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Semi_Direct) Call Close_SemiDSCF
      FstItr=.False.
!
      Call mma_deallocate(ip_ij)
      Call mma_deallocate(DMax)
      Call mma_deallocate(TMax)
!
      Call Term_Ints()
!
      Call Free_DeDe(Dens,TwoHam,nDens)
      Int_PostProcess => Null()
!
!
!     Broadcast contributions to the Fock matrix
!
      Call Sync_TH(TwoHam,nDens)
!                                                                      *
!***********************************************************************
!                                                                      *
!MAW start
!     CALL fmm_call_get_J_matrix(nDens,1,dens,TwoHam)
!MAW end
      Call Free_iSD()
!     Call Init_Int_Options()    ?
      End

      Subroutine Init_SemiDSCF(FstItr,Thize,Cutint)
      use dEAF
      use IOBUF, only: IODone, Disk, iBuf, ipos, iStatIO, Mode_Write,
     &                 OnDisk, Mode_Read, Disk_1, Disk_2, lBuf, nBuf,
     &                 Buffer, ID, LuTmp
      implicit None
#include "SysDef.fh"
      Logical FstItr
      Real*8 Thize, CutInt

      Integer lBufOld, nBufOld
      Real*8 ThizeOld, CutIntOld
      real*8 control(4)
!     Write (6,*) 'Enter: Init_SemiDSCF'
!     Write (6,*) 'Ondisk=',Ondisk
!     Write (6,*) 'lBuf=',lBuf
!
!---- Initiate asynchronous double buffer I/O.
!
      IODone = .False.
      Disk = 0.0D0
      iBuf=1
      iPos = 1
      If (FstItr) Then
         iStatIO = Mode_Write
!        write(6,*) 'write istatio=',istatio
         control(1)=Dble(lbuf)
         control(2)=Dble(nbuf)
         control(3)=thize
         control(4)=cutint
!        write(6,*) 'control written:',control
!        Write (6,*) ' Initiate write @', Disk,'iBuf=',iBuf
         If(OnDisk) Call dEAFAwrite(LuTmp,control,4*RtoI,Disk,id)
      Else
         iStatIO = Mode_Read
!        write(6,*) 'read istatio=',istatio
!
!------- Initiate first read ahead of time.
!
!        Write (6,*) 'lBuf*RtoI=',lbuf*RtoI,' rtoi=',Rtoi
         If (OnDisk) then
!           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
            Call dEAFread(LuTmp,control,4*RtoI,Disk)
            Disk_2 = Disk
            Disk_1 = Disk
!           write(6,*) 'control read:',control
            lbufold=nint(control(1))
            nbufold=nint(control(2))
            thizeold=control(3)
            cutintold=control(4)
            if (lbufold.lt.lbuf) then
              write(6,*) 'Reducing the buffer size from ',lbuf,
     &                  ' to',lbufold
              lbuf=lbufold
            else if(lbufold.gt.lbuf) then
              write(6,*) 'Inconsistent buffer lengths. Old:',lbufold,
     &                   '  current:',lbuf
              call Abend()
            end if
            if(nbuf.ne.nbufold) then
              write(6,*) 'Inconsistent buffer number. Old:',nbufold,
     &                   '  current:',nbuf
              call Abend()
            end if
            if(abs(thize-thizeold).gt.1.d-10) then
              write(6,*) 'Resetting thize from',thize,' to',thizeold
              thize=thizeold
            end if
            if(cutintold.gt.cutint) then
              write(6,*) 'Inconsistent Cutint. Old:',cutintold,
     &                   '  current:',cutint
              call Abend()
            end if
!           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
!           If(OnDisk) Write (6,*) ' Initial EAFARead'
            Call dEAFARead(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk,id)
         End If
      End If
!
!     Write (*,*) 'Exit: Init_SemiDSCF'
      End Subroutine Init_SemiDSCF

      Subroutine Close_SemiDSCF()
      use IOBUF, only: iPos, OnDisk, iStatIO, iBuf, lStRec, Mode_None
      Implicit None
!     Write (6,*) 'Enter: Close_SemiDSCF'
!
!---- If data was transfered to the I/O buffer write buffer on disc.
!
!  If buffer empty force the write :
      If (iPos.EQ.1) iPos=2
      If (OnDisk) Call WLBuf
!
      iPos = lStRec+1
      iStatIO = Mode_None
      iBuf = -99
!     Write (6,*) 'Exit: Close_SemiDSCF'
!
      End Subroutine Close_SemiDSCF

      Subroutine Mode_SemiDSCF(Wr_Mode)
      use IOBUF, only: iStatIO, Mode_Read, Disk, Disk_2, Mode_Write
      Implicit None
      Logical Wr_Mode
!
!     Write (6,*) 'Mode_SemiDSCF: Wr_Mode=',Wr_Mode
      If (Wr_Mode) Then
         If (iStatIO.eq.Mode_Read) Then
            Disk = Disk_2
            iStatIO = Mode_Write
!           Write (6,*) 'Changing to Write mode @',Disk
         End If
      Else
         If (iStatIO.eq.Mode_Write) Then
            Write (6,*) 'Change from Write to Read mode not implemented'
            Call Abend()
         End If
      End If
!
      End Subroutine Mode_SemiDSCF
