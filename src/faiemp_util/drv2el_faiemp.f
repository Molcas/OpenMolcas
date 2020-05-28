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
* Copyright (C) Ben Swerts                                             *
*               2016, Liviu Ungur                                      *
************************************************************************
      SubRoutine Drv2El_FAIEMP()
************************************************************************
*                                                                      *
*  Object: driver for the central-fragment interaction 2-electron      *
*          integrals (based on drv2el_3center_RI and drv2el_scf)       *
*                                                                      *
* Called from: Drv1El                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              Timing                                                  *
*              Setup_Ints                                              *
*              Eval_Ints                                               *
*              Term_Ints                                               *
*              QExit                                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*   Modified: Liviu Ungur                                              *
************************************************************************
      use k2_arrays, only: FT, nFT
      Implicit None
      External No_Routine
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "setup.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer      nTInt
      Parameter(   nTInt=1)
      Real*8       TInt(nTInt)
*
      Logical      W2Disc, PreSch, FreeK2, Verbose, Indexation,
     &             DoIntegrals, DoFock, DoGrad,NoCoul,NoExch
      Integer      iTOffs(8,8,8),
     &             nShi(8), nShj(8), nShk(8), nShl(8),
     &             nShOffi(8), nShOffj(8), nShOffk(8), nShOffl(8)
      Integer      nBas_Valence(0:7)
      Character*8  Label
      Character*80 Line
      Logical      lNoSkip, EnergyWeight
      Integer      i, j, iCnt, iCnttp, iDpos, iFD, iFpos, iIrrep, ijS,
     &             Ind, iOpt, ip_ij, ipDens, ipDMax, ipDq, ipFock, ipFq,
     &             ipFragDensAO, ipOneHam, ipTMax, iRC, iPrint, iRout,
     &             ipFragDensSO, iS, jS, lS, kS, klS, maxDens, mdc,
     &             lOper, mDens, nBasC, nBT, nBVT, nBVTi, nFock, nij,
     &             nInd, nOneHam, Nr_Dens, nSkal, nSkal_Fragments,
     &             nSkal_Valence
      Dimension    Ind(1,1,2)

      Real*8       Aint, Count, Disc, Disc_Mx, Dix_Mx, Dtst, ExFac,
     &             P_Eff, TCpu1, TCpu2, Thize, ThrAO, TMax_all,
     &             TskHi, TskLw, TWall1, TWall2, DMax, TMax
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*----- Statement functions
*
      TMax(i,j)=Work((j-1)*nSkal+i+ipTMax-1)
      DMax(i,j)=Work((j-1)*nSkal+i+ipDMax-1)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 203
      iPrint = nPrint(iRout)
      Call QEnter('Drv2ElFrag')
      call xFlush(6)
      nInd=1
      ExFac=One
      Nr_Dens=1
      DoIntegrals=.False.
      DoGrad=.False.
      NoCoul=.False.
      NoExch=.False.
c     W2Disc=.False.
      W2Disc=.True.
*
*     Handle both the valence and the fragment basis set
*
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal_Valence)
      Call Free_iSD
      Call Set_Basis_Mode('WithFragments')
      Call SetUp_iSD
      nBT = 0
      nBVT = 0
      do i = 0, nIrrep - 1
        nBas_Valence(i) = nBas(i)
        nBVT = nBVT + nBas(i)*(nBas(i)+1)/2
        nBas(i) = nBas(i) + nBas_Frag(i)
        nBT = nBT + nBas(i)*(nBas(i)+1)/2
c       write(*,*) 'For irrep ',i,', nBas_Valence, nBas_Frag, nBas =',
c     &                            nBas_Valence(i), nBas_Frag(i), nBas(i)
      enddo
c     write(*,*) 'nBT, nBVT =',nBT, nBVT
c     write(*,*) 'iOper =',(iOper(i),i=0,nIrrep-1)
c     mdc = 0
c     Do iCnttp = 1, nCnttp
c       Do iCnt = 1, nCntr(iCnttp)
c         mdc = mdc + 1
c         write(*,*) 'For iCnttp, iCnt, mdc:',iCnttp,iCnt,mdc
c         write(*,*) '   iChCnt  =',iChCnt(mdc),'(',
c    &     (iAnd(iChCnt(mdc),i),i=0,nIrrep-1),')'
c         write(*,*) '   jStab() =',(jStab(i,mdc),i=0,nIrrep-1)
c        write(*,*) 'nStab, nIrrep/nStab =',nStab(mdc),nIrrep/nStab(mdc)
c         Do iIrrep = 0, nIrrep - 1
c           If(iAnd(iChCnt(mdc),iIrrep).eq.iOper(iIrrep))
c    &        write(*,*) '               center appears in irrep',iIrrep
c         End Do
c       End Do
c     End Do
*                                                                      *
************************************************************************
*                                                                      *
*---  Construct custom density matrix
*

      Call GetMem('Density','Allo','Real',ipDens,nBT)
      Call GetMem('Result','Allo','Real',ipFock,nBT)
* Valence part is zero
      call dcopy_(nBT, [Zero], 0, Work(ipDens), 1)
      call dcopy_(nBT, [Zero], 0, Work(ipFock), 1)
* Each fragment needs it's (symmetrized) density matrix added along the diagonal
* This density matrix first has to be constructed from the MO coefficients
* so allocate space for the largest possible density matrix
      maxDens = 0
      Do iCnttp = 1, nCnttp
        If(nFragType(iCnttp).gt.0) maxDens = Max(maxDens,
     &                        nFragDens(iCnttp)*(nFragDens(iCnttp)+1)/2)
      End Do
c      write(*,*) 'maxDens =',maxDens
      Call GetMem('FragDSO','Allo','Real',ipFragDensSO,maxDens)
c      If(nIrrep.ne.1) Then
c        Call GetMem('FragDAO','Allo','Real',ipFragDensAO,maxDens)
c      Else
        ipFragDensAO = ipFragDensSO
c      End If

      iDpos = 1 ! position in the total density matrix
      Do iIrrep = 0, nIrrep - 1
        nBasC = nBas_Valence(iIrrep)
        iDpos = iDpos + nBasC*(nBasC+1)/2
        mdc = 0
        Do 1000 iCnttp = 1, nCnttp
          If(nFragType(iCnttp).le.0) Then
            mdc = mdc + nCntr(iCnttp)
            Go To 1000
          End If
* construct the density matrix
          EnergyWeight = .false.
c          write(6,*) 'Drv2ElFrag: call MakeDens:',nFragDens(iCnttp),
c     &                                            nFragEner(iCnttp)
c          call xFlush(6)
c         Call MakeDens(nFragDens(iCurCnttp),nFragEner(iCurCnttp),
c     &                Work(ipFragCoef(iCurCnttp)),Work(ipFragEner(iCurCnttp)),
c     &                EnergyWeight,Array)
          Call MakeDens(nFragDens(iCnttp),nFragEner(iCnttp),
     &                Work(ipFragCoef(iCnttp)),Work(ipFragEner(iCnttp)),
     &                EnergyWeight,Work(ipFragDensAO))
* create the symmetry adapted version if necessary
* (fragment densities are always calculated without symmetry)
C         If(nIrrep.ne.1) Call SymmDens(Work(ipFragDensAO),
C    &      Work(ipFragDensSO))
          If(iPrint.ge.99) Call TriPrt('Fragment density',' ',
     &      Work(ipFragDensSO),nFragDens(iCnttp))

          Do iCnt = 1, nCntr(iCnttp)
            mdc = mdc + 1
* only add fragment densities that are active in this irrep
* => the following procedure still has to be verified thoroughly
*    but appears to be working
            If(iAnd(iChCnt(mdc),iIrrep).eq.iOper(iIrrep)) Then
* add it at the correct location in the large custom density matrix
              iFpos = 1
c              ! position in fragment density matrix
              Do i = 1, nFragDens(iCnttp)
                iDpos = iDpos + nBasC
                Do j = 0, i-1
                  Work(ipDens + iDpos + j - 1) =
     &                          Work(ipFragDensSO + iFpos + j - 1)
                End Do
                iDpos = iDpos + i
                iFpos = iFpos + i
              End Do
              nBasC = nBasC + nFragDens(iCnttp)
            End If
          End Do
 1000   Continue
      End Do
      If(iPrint.ge.19) Then
        iFD = ipDens
        Do iIrrep = 0, nIrrep - 1
          Call TriPrt('Combined density',' ',Work(iFD),nBas(iIrrep))
          iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
        End Do
      End If
      Call GetMem('FragDSO','Free','Real',ipFragDensSO,maxDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate auxiliary memory needed for direct integral evaluation
*     if the Fock matrix
*
      DoFock=.True.
      If (DoFock) Then
         nFT = MxFT
         Call mma_allocate(FT,nFT,Label='FT')
         MxDij = 6 * nIrrep * MxDij
         Call GetMem('Dijs','Allo','Real',ipDijs,MxDij)
      Else
         ipDijs=ip_iDummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Desymmetrize the custom density matrix
*
* Should be possible to reduce the storage space for the Fock matrix
* and only store the top nBas_Valence(0) rows (non-symmetry case tested)
*
c      write(6,*) 'Drv2ElFrag: call AlloK2:'
c      call xFlush(6)
      Call AlloK2()
c      write(6,*) 'Drv2ElFrag: call DEDE2_FAIEMP: nBT,mDens: ',nBT,mDens
c      call xFlush(6)
      Call DeDe_FAIEMP(Work(ipDens),Work(ipFock),nBT,mDens,ipDq,ipFq)
c     write(*,*) 'Drv2El_FAIEMP: mDens =',mDens
      If(iPrint.ge.99) Then
        If(nIrrep.eq.1) Then
          Call RecPrt('Desymmetrized Density:',' ',Work(ipDq),nBas(0),
     &                nBas(0))
        Else
          iFD = ipDens
          Do iIrrep = 0, nIrrep - 1
         Call TriPrt('Desymmetrized density',' ',Work(iFD),nBas(iIrrep))
            iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
          End Do
        End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      ThrAO = Zero
      Indexation = .False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      nSkal_Fragments=nSkal-nSkal_Valence
c      write(*,*) 'Drv2El_FAIEMP: nSkal =',nSkal,' (nSkal_Valence =',
c     &        nSkal_Valence,' and nSkal_Fragments =',nSkal_Fragments,')'
*                                                                      *
************************************************************************
*                                                                      *
      Thize=1.0d-6
      PreSch=.False.
      Disc_Mx=Zero
*
      Disc=Zero
      Dix_Mx=Zero
      TskHi=Zero
      TskLw=Zero
      ThrInt = CutInt   ! Integral neglect threshold from SCF
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
c     write(*,*) 'Drv2El_FAIEMP: TMax_All =',TMax_All
      Call GetMem('DMax','Allo','Real',ipDMax,nSkal**2)
      Call Shell_MxDens(Work(ipDens),work(ipDMax),nSkal)
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
      P_Eff=dble(nij)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*
*     Now do a quadruple loop over shells
*
c     ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
      ijS = 1
      iS = iWork((ijS-1)*2+ip_ij)
      jS = iWork((ijS-1)*2+ip_ij+1)
c     klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      klS = 1
      kS = iWork((klS-1)*2+ip_ij)
      lS = iWork((klS-1)*2+ip_ij+1)
 13   Continue
      If(ijS.gt.int(P_Eff)) Go To 12
*                                                                      *
************************************************************************
*                                                                      *
* density prescreening (results in iS > nSkal_Valence)
         Aint=TMax(iS,jS)*TMax(kS,lS)
         Dtst=Max(DMax(is,ls)/Four,DMax(is,ks)/Four,
     &            DMax(js,ls)/Four,DMax(js,ks)/Four,
     &            DMax(is,js),DMax(ks,ls))
         lNoSkip = Aint*Dtst.ge.ThrInt
* only calculate needed integrals and only update the valence part of the
* Fock matrix (iS > nSkal_Valence, lS <= nSkal_Valence, jS and kS
* belonging to different regions)
         If(jS.le.nSkal_Valence) Then
           lNoSkip = lNoSkip.and.kS.gt.nSkal_Valence
         Else
           lNoSkip = lNoSkip.and.kS.le.nSkal_Valence
         End If
         lNoSkip = lNoSkip.and.lS.le.nSkal_Valence

         If (lNoSkip) Then
           Call Eval_Ints_New_
     &                    (iS,jS,kS,lS,TInt,nTInt,
     &                     iTOffs,nShi,nShj,nShk,nShl,
     &                     nShOffi,nShOffj,nShOffk,nShOffl,
     &                     No_Routine,
     &                     Work(ipDq),Work(ipFq),mDens,[ExFac],Nr_Dens,
     &                     Ind,nInd,[NoCoul],[NoExch],
     &                     Thize,W2Disc,PreSch,Disc_Mx,Disc,
     &                     Count,DoIntegrals,DoFock)
           If(iPrint.ge.99) Then
            write(6,*) 'Drv2El_FAIEMP: for iS, jS, kS, lS =',is,js,ks,ls
             If(nIrrep.eq.1) Then
              Call RecPrt('updated Fock',' ',Work(ipFq),nBas(0),nBas(0))
             Else
               iFD = ipFq
               Do iIrrep = 0, nIrrep - 1
                 Call TriPrt('updated Fock',' ',Work(iFD),nBas(iIrrep))
                 iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
               End Do
             End If
           End If
         End If
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
*     Use a time slot to save the number of tasks and shell
*     quadruplets processed by an individual node
      Call SavStat(1,One,'+')
      Call SavStat(2,TskHi-TskLw+One,'+')
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('ip_ij','Free','Inte',ip_ij,nSkal*(nSkal+1))
      Call GetMem('DMax','Free','Real',ipDMax,nSkal**2)
      Call GetMem('TMax','Free','Real',ipTMax,nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*
      Call Free_DeDe_FAIEMP(Work(ipDens),Work(ipFock),nBT,ipDq,ipFq)

      Call GetMem('Density','Free','Real',ipDens,nBT)
      If(iPrint.ge.10) Then
        write(6,*)
        write(6,*)
        write(6,'(a)') 'SO Integrals of type Frag2El Component 1'
        iFD = ipFock
        Do iIrrep = 0, nIrrep - 1
          Write (Line,'(1X,A,I1)')
     &      ' Diagonal Symmetry Block ', iIrrep+1
          Call TriPrt(Line,' ',Work(iFD),nBas_Valence(iIrrep))
          iFD = iFD + nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
        End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write the results to the one electron integral file
*
* read the one electron hamiltonian
      Label = 'OneHam  '
      iRC = -1
      iOpt = 0
      Call GetMem('Temp','Allo','Real',ipOneHam,nBVT+4)
      Call RdOne(iRC,iOpt,Label,1,Work(ipOneHam),lOper)
      If(iRC.ne.0) Then
        Write (6,*) 'Drv2El_FAIEMP: Error reading from ONEINT'
        Write (6,'(A,A)') 'Label=',Label
        Call Abend()
      End If
c     If(nIrrep.eq.1) Then
c       Call TriPrt('OneHam as read from OneInt',' ',Work(ipOneHam),
c    &              nBas_Valence(0))
c     Else
c       iFD = ipOneHam
c       Do iIrrep = 0, nIrrep - 1
c        Call TriPrt('OneHam before',' ',Work(iFD),nBas_Valence(iIrrep))
c         iFD = iFD + nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
c       End Do
c     End If
* add the calculated results
      nOneHam = 0 ! counter in the ipOneHam matrices (small)
      nFock = 0   ! counter in the ipFock matrices (larger)
      Do iIrrep = 0, nIrrep - 1
        nBVTi = nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
        call daxpy_(nBVTi,One,Work(ipFock+nFock),1,
     &                       Work(ipOneHam+nOneHam),1)
        nOneHam = nOneHam + nBVTi
        nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do

* write out the results
      If(iPrint.ge.15) Then
        iFD = ipOneHam
        Do iIrrep = 0, nIrrep - 1
        Call TriPrt('OneHam at end',' ',Work(iFD),nBas_Valence(iIrrep))
          iFD = iFD + nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
        End Do
      End If
      iRC = -1
      Call WrOne(iRC,iOpt,Label,1,Work(ipOneHam),lOper)
      If(iRC.ne.0) Then
        Write (6,*) 'Drv2El_FAIEMP: Error writing to ONEINT'
        Write (6,'(A,A)') 'Label=',Label
        Call Abend()
      End If
      iRC = -1
      Label = 'OneHam 0'
      Call WrOne(iRC,iOpt,Label,1,Work(ipOneHam),lOper)
      If(iRC.ne.0) Then
        Write (6,*) 'Drv2El_FAIEMP: Error writing to ONEINT'
        Write (6,'(A,A)') 'Label=',Label
        Call Abend()
      End If

* cleanup
      Call GetMem('Temp','Free','Real',ipOneHam,nBVT+4)
      Call GetMem('Result','Free','Real',ipFock,nBT)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      Call Free_iSD
      Call Set_Basis_Mode('Valence')
      Call SetUp_iSD
      do i = 0, nIrrep - 1
      nBas(i) = nBas_Valence(i)
      enddo
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Drv2ElFrag')
c      write(6,*) 'Exiting Drv2ElFrag'
c      call xFlush(6)

      Return
      End
