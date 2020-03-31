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
* Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
*               1990, IBM                                              *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine Drv2El_dscf(Dens,TwoHam,nDens,nDisc,Thize,PreSch,
     &                       FstItr,NoCoul,ExFac)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will control the type of the two-electron integral, e.g.    *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis func-  *
*          tions of the requested type.                                *
*                                                                      *
*          Dens is the folded lower triangular of the 1st order        *
*               density matrix.                                        *
*          Twoham is the lower triangular of the two-electron contri-  *
*               bution to the Fock matrix.                             *
*                                                                      *
* Called from: PMat                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              DeDe_SCF                                                *
*              DrvK2                                                   *
*              StatP                                                   *
*              GetMem                                                  *
*              mHrr                                                    *
*              DCopy   (ESSL)                                          *
*              Swap                                                    *
*              DrvTwo                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified for direct SCF. January '93                     *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             Modified by M.Schuetz @teokem.lu.se :                    *
*             parallel region split off in drvtwo.f, April '95         *
*             Modified by R. Lindh  @teokem.lu.se :                    *
*             total repacking of code September '96                    *
************************************************************************

      Implicit Real*8 (a-h,o-z)
      External Rsv_GTList, No_Routine
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "real.fh"
#include "k2.fh"
#include "nsd.fh"
#include "setup.fh"
      Logical NoCoul,NoExch
#include "IOBuf.fh"
*
      Parameter(nTInt=1)
      Real*8 Dens(nDens), TwoHam(nDens), TInt(nTInt)
      Logical W2Disc, FstItr, Semi_Direct,Rsv_GTList,
     &        PreSch, Density, Free_K2, Verbose, Indexation,
     &        DoIntegrals, DoFock, DoGrad, Triangular
      Integer iTOffs(8,8,8),
     &        nShi(8), nShj(8), nShk(8), nShl(8),
     &        nShOffi(8), nShOffj(8), nShOffk(8), nShOffl(8)
      Logical Debug
      Character*72 SLine
      Dimension Ind(1,1,2)
*                                                                      *
************************************************************************
*                                                                      *
*---- Statement functions
*
      TMax(i,j)=Work((j-1)*nSkal+i+ipTMax-1)
      DMax(i,j)=Work((j-1)*nSkal+i+ipDMax-1)

*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('Drv2El')
#ifdef _DEBUG_
       Debug=.true.
c       iPrint=200
#else
       Debug=.false.
#endif
*
      SLine='Computing 2-electron integrals'
      Call StatusLine(' SCF:',SLine)
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      nInd=1
      Nr_Dens=1
      DoIntegrals=.False.
      DoFock=.True.
      DoGrad=.False.
      NoExch=ExFac.eq.Zero
*                                                                      *
************************************************************************
*                                                                      *
*-----Set up for partial SO/AO integral storage.
*
*---- nDisc = file size in kbyte from input
      Semi_Direct = nDisc.ne.0
      If (Semi_Direct) Call Init_SemiDSCF(FstItr,Thize,Cutint)
*     Disc_Mx = file size in Real*8 128=1024/8
      Disc_Mx= DBLE(nDisc)*128.D00
*     Subtract for the last buffer
      Disc_Mx= Disc_Mx - lBuf
*                                                                      *
************************************************************************
*                                                                      *
*-----Desymmetrize differential densities.
*     Observe that the desymmetrized 1st order density matrices are
*     canonical, i.e. the relative order of the indices are canonically
*     ordered.
*
      Density=.True.       ! Use density information in prescreening.
      Call DeDe_SCF(Dens,TwoHam,nDens,mDens,ipDq,ipFq)

*                                                                      *
************************************************************************
*                                                                      *
      ThrAO=Zero           ! Do not modify CutInt
      Indexation=.False.
*
      Call SetUp_Ints(nSkal,Indexation,ThrAO,Density,DoGrad)
*                                                                      *
************************************************************************
*                                                                      *
      Disc = Zero
      W2Disc=.False.
      TskHi=Zero
      TskLw=Zero
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call GetMem('TMax','Allo','Real',ipTMax,nSkal**2)
      Call Shell_MxSchwz(nSkal,Work(ipTMax))
      TMax_all=Zero
      Do iS = 1, nSkal
         Do jS = 1, iS
            TMax_all=Max(TMax_all,TMax(iS,jS))
         End Do
      End Do
      Call GetMem('DMax','Allo','Real',ipDMax,nSkal**2)
      Call Shell_MxDens(Dens,work(ipDMax),nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      Call GetMem('ip_ij','Allo','Inte',ip_ij,nSkal*(nSkal+1))
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nij = nij + 1
               iWork((nij-1)*2+ip_ij  )=iS
               iWork((nij-1)*2+ip_ij+1)=jS
            End If
         End Do
      End Do
      P_Eff=Dble(nij)
*
      PP_Eff=P_Eff**2
      PP_Eff_delta=0.10D0*PP_Eff
      PP_Count=Zero
*                                                                      *
************************************************************************
*                                                                      *
*.... For distributed parallel SCF initiate (sequential code is special
*     case when the number of nodes in the mpp is 1).
*
*     1: Task list (tlist)
*     2: Private priority list (pplist)
*     3: Global task list (gtlist)
*
      If (FstItr) Then
         Triangular=.True.
         Call Alloc_TList(Triangular,P_Eff)
         Call Init_TList(Triangular,P_Eff)
         Call Init_PPList
         Call Init_GTList
      Else
         Call ReInit_PPList(Semi_Direct)
         Call ReInit_GTList
      End If
      iOpt=0
      If (.Not.FstItr.and.Semi_direct) iOpt=2
*
      Call CWTime(TCpu1,TWall1)
*
*     big loop over individual tasks, distributed over individual nodes

   10 Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_GTList(TskLw,TskHi,iOpt,W2Disc)) Go To 11
      Call Mode_SemiDSCF(W2Disc)
*     Write (6,*) 'TskLw,TskHi,W2Disc=',TskLw,TskHi,W2Disc
*
*     Now do a quadruple loop over shells
*
      ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
      iS = iWork((ijS-1)*2+ip_ij)
      jS = iWork((ijS-1)*2+ip_ij+1)
      klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      kS = iWork((klS-1)*2+ip_ij)
      lS = iWork((klS-1)*2+ip_ij+1)
      Count=TskLw
      If (Count-TskHi.gt.1.0D-10) Go To 12
  13  Continue
*
      S_Eff=DBLE(ijS)
      T_Eff=DBLE(klS)
      ST_Eff=S_Eff*(S_Eff-One)/2D0 + T_Eff
      If (ST_Eff.ge.PP_Count) Then
         Write (SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',
     &        ST_Eff/PP_Eff,'% done so far.'
         Call StatusLine(' Seward:',SLine)
         PP_Count = PP_Count + PP_Eff_delta
      End If
*                                                                      *
************************************************************************
*                                                                      *
         Aint=TMax(iS,jS)*TMax(kS,lS)
         If (Semi_Direct) Then
*
*           No density screening in semi-direct case!
*           Cutint: Threshold for Integrals. In semi-direct case, this
*                   must be the final threshold used in the last scf
*                   iteration
*           Thrint: Threshold for Density screening. This the actual
*                   threshold
*                   for the current iteration
*
           If (AInt.lt.CutInt) Go To 14
         Else
           If(NoCoul) then
              Dtst=Max(DMax(is,ls)/Four,DMax(is,ks)/Four,
     &                 DMax(js,ls)/Four,DMax(js,ks)/Four)
           Else If(NoExch) then
              Dtst=Max(DMax(is,js),DMax(ks,ls))
           Else
              Dtst=Max(DMax(is,ls)/Four,DMax(is,ks)/Four,
     &                 DMax(js,ls)/Four,DMax(js,ks)/Four,
     &                 DMax(is,js),DMax(ks,ls))
           End If
           If (Aint*Dtst.lt.ThrInt) goto 14
         End if
*                                                                      *
************************************************************************
*                                                                      *
         Call Eval_Ints_New_
     &                  (iS,jS,kS,lS,TInt,nTInt,
     &                   iTOffs,nShi,nShj,nShk,nShl,
     &                   nShOffi,nShOffj,nShOffk,nShOffl,
     &                   No_Routine,
     &                   Work(ipDq),Work(ipFq),mDens,[ExFac],Nr_Dens,
     &                   Ind,nInd,[NoCoul],[NoExch],
     &                   Thize,W2Disc,PreSch,Disc_Mx,Disc,
     &                   Count,DoIntegrals,DoFock)

*
 14      Continue
         Count=Count+One
         If (Count-TskHi.gt.1.0D-10) Go To 12
         klS = klS + 1
         If (klS.gt.ijS) Then
            ijS = ijS + 1
            klS = 1
         End If
         iS = iWork((ijS-1)*2+ip_ij  )
         jS = iWork((ijS-1)*2+ip_ij+1)
         kS = iWork((klS-1)*2+ip_ij  )
         lS = iWork((klS-1)*2+ip_ij+1)
         Go To 13
*
*     Task endpoint
*
 12   Continue
*
      If (Semi_Direct) Then
         If (W2Disc) Then
            Call Put_QLast
         Else
            Call Pos_QLast(Disc)
         End If
      End If
*
*     Use a time slot to save the number of tasks and shell
*     quadruplets processed by an individual node
      Call SavStat(1,One,'+')
      Call SavStat(2,TskHi-TskLw+One,'+')
      Go To 10
 11   Continue
*     End of big task loop
      Lu=6
      Call CWTime(TCpu2,TWall2)
      Call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      If (Semi_Direct) Call Close_SemiDSCF
      FstItr=.False.
*
      Call GetMem('ip_ij','Free','Inte',ip_ij,nSkal*(nSkal+1))
      Call GetMem('DMax','Free','Real',ipDMax,nSkal**2)
      Call GetMem('TMax','Free','Real',ipTMax,nSkal**2)
*
      Verbose=.False.
      Free_K2=.False. ! Call to freek2 is external to the driver.
      Call Term_Ints(Verbose,Free_K2)
*
      Call Free_DeDe2(Dens,TwoHam,nDens,ipDq,ipFq)
*
      Call QExit('Drv2El')
*
*     Broadcast contributions to the Fock matrix
*
      Call Sync_TH(TwoHam,nDens)
*                                                                      *
************************************************************************
*                                                                      *
CMAW start
#ifdef _F90ENABLE_
      CALL fmm_call_get_J_matrix(nDens,1,dens,TwoHam)
#endif
CMAW end
      Call Free_iSD()
      Return
      End
      Subroutine Init_SemiDSCF(FstItr,Thize,Cutint)
      use dEAF
      implicit real*8 (a-h,o-z)
#include "IOBuf.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
      real*8 control(4)
      Logical FstItr
*     Write (6,*) 'Enter: Init_SemiDSCF'
*     Write (6,*) 'Ondisk=',Ondisk
*     Write (6,*) 'lBuf=',lBuf
*
*---- Initiate asynchronous double buffer I/O.
*
      IODone = .False.
      Disk = 0.0D0
      iBuf=1
      iPos = 1
      If (FstItr) Then
         iStatIO = Mode_Write
*        write(6,*) 'write istatio=',istatio
         control(1)=Dble(lbuf)
         control(2)=Dble(nbuf)
         control(3)=thize
         control(4)=cutint
*        write(6,*) 'control written:',control
C        Write (6,*) ' Initiate write @', Disk,'iBuf=',iBuf
         If(OnDisk) Call dEAFAwrite(LuTmp,control,4*RtoI,Disk,id)
      Else
         iStatIO = Mode_Read
*        write(6,*) 'read istatio=',istatio
*
*------- Initiate first read ahead of time.
*
*        Write (6,*) 'lBuf*RtoI=',lbuf*RtoI,' rtoi=',Rtoi
*        Call GetMem('ISemi','List','Real',iDum,iDum)
         If (OnDisk) then
C           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
            Call dEAFread(LuTmp,control,4*RtoI,Disk)
            Disk_2 = Disk
            Disk_1 = Disk
*           write(6,*) 'control read:',control
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
              call Abend
            end if
            if(nbuf.ne.nbufold) then
              write(6,*) 'Inconsistent buffer number. Old:',nbufold,
     &                   '  current:',nbuf
              call Abend
            end if
            if(abs(thize-thizeold).gt.1.d-10) then
              write(6,*) 'Resetting thize from',thize,' to',thizeold
              thize=thizeold
            end if
            if(cutintold.gt.cutint) then
              write(6,*) 'Inconsistent Cutint. Old:',cutintold,
     &                   '  current:',cutint
              call Abend
            end if
c           Write (6,*) ' Initiate read @', Disk,'iBuf=',iBuf
*           If(OnDisk) Write (6,*) ' Initial EAFARead'
            Call dEAFARead(LuTmp,Work((iBuf-1)*lBuf+ipBuf),
     &                               lBuf*RtoI,Disk,id)
         End If
      End If
*
*     Write (*,*) 'Exit: Init_SemiDSCF'
      Return
      End
      Subroutine Close_SemiDSCF
#include "IOBuf.fh"
#include "WrkSpc.fh"
*     Write (6,*) 'Enter: Close_SemiDSCF'
*
*---- If data was transfered to the I/O buffer write buffer on disc.
*
*     Call GetMem('CSemi','List','Real',iDum,iDum)
C  If buffer empty force the write :
      If (iPos.EQ.1) iPos=2
      If (OnDisk) Call WLBuf
*
      iPos = lStRec+1
      iStatIO = Mode_None
      iBuf = -99
*     Write (6,*) 'Exit: Close_SemiDSCF'
*
      Return
      End
      Subroutine Mode_SemiDSCF(Wr_Mode)
#include "IOBuf.fh"
      Logical Wr_Mode
*
*     Write (6,*) 'Mode_SemiDSCF: Wr_Mode=',Wr_Mode
      If (Wr_Mode) Then
         If (iStatIO.eq.Mode_Read) Then
            Disk = Disk_2
            iStatIO = Mode_Write
*           Write (6,*) 'Changing to Write mode @',Disk
         End If
      Else
         If (iStatIO.eq.Mode_Write) Then
            Write (6,*) 'Change from Write to Read mode not implemented'
            Call Abend
         End If
      End If
*
      Return
      End
