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
      use k2_arrays, only: pDq, pFq
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      use Real_Info, only: ThrInt, CutInt
      Implicit None
      External No_Routine
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
      Integer      iTOffs(8,8,8)
      Integer      nBas_Valence(0:7)
      Character*8  Label
      Logical      lNoSkip, EnergyWeight
      Integer      i, j, iCnt, iCnttp, iDpos, iFpos, iIrrep, ijS,
     &             iOpt, ip_ij, ipDMax,
     &             ipFragDensAO, ipOneHam, ipTMax, iRC, iPrint, iRout,
     &             ipFragDensSO, iS, jS, lS, kS, klS, maxDens, mdc,
     &             lOper, mDens, nBasC, nBT, nBVT, nBVTi, nFock, nij,
     &             nOneHam, Nr_Dens, nSkal, nSkal_Fragments,
     &             nSkal_Valence

      Real*8       Aint, Count, Disc, Disc_Mx, Dix_Mx, Dtst, ExFac,
     &             P_Eff, TCpu1, TCpu2, Thize, ThrAO, TMax_all,
     &             TskHi, TskLw, TWall1, TWall2, DMax, TMax
      Real*8, Allocatable, Target:: Dens(:), Fock(:)
*define _DEBUG_
#ifdef _DEBUG_
      Integer      iFD
      Character*80 Line
#endif
#include "../integral_util/dede_interface.fh"
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
      call xFlush(6)
      ExFac=One
      Nr_Dens=1
      DoIntegrals=.False.
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
      enddo
*                                                                      *
************************************************************************
*                                                                      *
*---  Construct custom density matrix
*
      Call mma_allocate(Dens,nBT,Label='Dens')
      Call mma_allocate(Fock,nBT,Label='Fock')
* Valence part is zero
      Dens(:)=Zero
      Fock(:)=Zero
* Each fragment needs it's (symmetrized) density matrix added along the
* diagonal.
* This density matrix first has to be constructed from the MO coeffs
* so allocate space for the largest possible density matrix
      maxDens = 0
      Do iCnttp = 1, nCnttp
        If(dbsc(iCnttp)%nFragType.gt.0) maxDens = Max(maxDens,
     &     dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
      End Do
      Call GetMem('FragDSO','Allo','Real',ipFragDensSO,maxDens)
      ipFragDensAO = ipFragDensSO

      iDpos = 1 ! position in the total density matrix
      Do iIrrep = 0, nIrrep - 1
        nBasC = nBas_Valence(iIrrep)
        iDpos = iDpos + nBasC*(nBasC+1)/2
        mdc = 0
        Do 1000 iCnttp = 1, nCnttp
          If(dbsc(iCnttp)%nFragType.le.0) Then
            mdc = mdc + dbsc(iCnttp)%nCntr
            Go To 1000
          End If
* construct the density matrix
          EnergyWeight = .false.
          Call MakeDens(dbsc(iCnttp)%nFragDens,dbsc(iCnttp)%nFragEner,
     &                dbsc(iCnttp)%FragCoef,dbsc(iCnttp)%FragEner,
     &                EnergyWeight,Work(ipFragDensAO))
* create the symmetry adapted version if necessary
* (fragment densities are always calculated without symmetry)
#ifdef _DEBUG_
          Call TriPrt('Fragment density',' ',
     &      Work(ipFragDensSO),dbsc(iCnttp)%nFragDens)
#endif

          Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
* only add fragment densities that are active in this irrep
* => the following procedure still has to be verified thoroughly
*    but appears to be working
            If(iAnd(dc(mdc)%iChCnt,iIrrep).eq.iOper(iIrrep)) Then
* add it at the correct location in the large custom density matrix
              iFpos = 1
c              ! position in fragment density matrix
              Do i = 1, dbsc(iCnttp)%nFragDens
                iDpos = iDpos + nBasC
                Do j = 0, i-1
                  Dens(iDpos + j) =
     &                          Work(ipFragDensSO + iFpos + j - 1)
                End Do
                iDpos = iDpos + i
                iFpos = iFpos + i
              End Do
              nBasC = nBasC + dbsc(iCnttp)%nFragDens
            End If
          End Do
 1000   Continue
      End Do
#ifdef _DEBUG_
      iFD = 1
      Do iIrrep = 0, nIrrep - 1
         Call TriPrt('Combined density',' ',Dens(iFD),nBas(iIrrep))
         iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
#endif
      Call GetMem('FragDSO','Free','Real',ipFragDensSO,maxDens)
*                                                                      *
************************************************************************
*                                                                      *
*-----Desymmetrize the custom density matrix
*
* Should be possible to reduce the storage space for the Fock matrix
* and only store the top nBas_Valence(0) rows (non-symmetry case tested)
*
      Call AlloK2()
      Call DeDe_SCF(Dens,Fock,nBT,mDens)
#ifdef _DEBUG_
      If (nIrrep.eq.1) Then
         Call RecPrt('Desymmetrized Density:',' ',pDq,nBas(0),nBas(0))
      Else
         iFD = 1
         Do iIrrep = 0, nIrrep - 1
            Call TriPrt('Desymmetrized density',' ',
     &                  pDq(iFD),nBas(iIrrep))
            iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      ThrAO = Zero
      Indexation = .False.
      DoFock=.True.
      DoGrad=.False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      nSkal_Fragments=nSkal-nSkal_Valence
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
           Call Eval_Ints_New_Internal
     &                    (iS,jS,kS,lS,TInt,nTInt,
     &                     iTOffs,No_Routine,
     &                     pDq,pFq,mDens,[ExFac],Nr_Dens,
     &                     [NoCoul],[NoExch],
     &                     Thize,W2Disc,PreSch,Disc_Mx,Disc,
     &                     Count,DoIntegrals,DoFock)
#ifdef _DEBUG_
            write(6,*) 'Drv2El_FAIEMP: for iS, jS, kS, lS =',is,js,ks,ls
            If (nIrrep.eq.1) Then
               Call RecPrt('updated Fock',' ',pFq,nBas(0),nBas(0))
            Else
               iFD = 1
               Do iIrrep = 0, nIrrep - 1
                 Call TriPrt('updated Fock',' ',pFq(iFD),nBas(iIrrep))
                 iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
               End Do
            End If
#endif
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
      Call Free_DeDe(Dens,Fock,nBT)

      Call mma_deallocate(Dens)
#ifdef _DEBUG_
      write(6,*)
      write(6,*)
      write(6,'(a)') 'SO Integrals of type Frag2El Component 1'
      iFD = 1
      Do iIrrep = 0, nIrrep - 1
         Write (Line,'(1X,A,I1)')
     &      ' Diagonal Symmetry Block ', iIrrep+1
         Call TriPrt(Line,' ',Fock(iFD),nBas_Valence(iIrrep))
         iFD = iFD + nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
      End Do
#endif
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
* add the calculated results
      nOneHam = 0 ! counter in the ipOneHam matrices (small)
      nFock = 1   ! counter in the ipFock matrices (larger)
      Do iIrrep = 0, nIrrep - 1
        nBVTi = nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
        call daxpy_(nBVTi,One,Fock(nFock),1,
     &                       Work(ipOneHam+nOneHam),1)
        nOneHam = nOneHam + nBVTi
        nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do

* write out the results
#ifdef _DEBUG_
      iFD = ipOneHam
      Do iIrrep = 0, nIrrep - 1
         Call TriPrt('OneHam at end',' ',Work(iFD),nBas_Valence(iIrrep))
         iFD = iFD + nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
      End Do
#endif
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
      Call mma_deallocate(Fock)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      Call Free_iSD
      Call Set_Basis_Mode('Valence')
      Call SetUp_iSD
      Do i = 0, nIrrep - 1
         nBas(i) = nBas_Valence(i)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
