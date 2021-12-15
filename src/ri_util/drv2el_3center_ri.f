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
* Copyright (C) 1990,1991,1993,1998,2006,2007, Roland Lindh            *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv2El_3Center_RI(Integral_WrOut,ThrAO)
************************************************************************
*                                                                      *
*  Object: driver for the 3 center integrals in the RI approach.       *
*                                                                      *
*          This code has three sections                                *
*          1) a 2-center section to generate the Q-vectors             *
*          2) a 3-center section to generate the R-vectors             *
*          3) a partial transpose section to generate the RI vectors   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             Modified driver. Jan. '98                                *
*             Modified to 3-center ERIs for RI Jan '06                 *
*             Modified to out-of-core version Feb '07                  *
************************************************************************
      use iSD_data
      use Wrj12
      use Basis_Info, only: dbsc, nBas, nBas_Aux
      use Temporary_Parameters, only: force_out_of_core
      use Real_Info, only: CutInt
      use RICD_Info, only: LDF
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut, Integral_RI_2, Rsv_Tsk
#include "Molcas.fh"
#include "j12.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
*
#include "lRI.fh"
#include "setup.fh"
#include "nsd.fh"
#define _no_nShs_
#include "iTOffs.fh"
*
#include "para_info.fh"
      Character*6 Name_R*6
      Integer iOff_3C(3,0:7), Lu_R(0:7), iAddr_R(0:7), iMax_R(2,0:7),
     &        iTtmp(0:7), NoChoVec(0:7), iOff_Rv(0:7)
      Logical Verbose, Indexation, FreeK2, Distribute,
     &        DoGrad, DoFock, Out_of_Core, Rsv_Tsk, Reduce_Prt
      External Reduce_Prt
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*----- Statement functions
*
      TMax_Valence(i,j)=Work(ipTMax-1+(j-1)*nSkal_Valence+i)
      TMax_Auxiliary(i)=Work(ipTMax-1+nSkal_Valence**2+i)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
*
#ifdef  _MOLCAS_MPP_
      Distribute = nProcs.gt.1 .and. Is_Real_Par()
#else
      Distribute = .False.
#endif
*
*     Get global print level
*
      iPL=iPrintLevel(-1)
      If (iPL.eq.2) Then
         iPL=5
      Else If (iPL.eq.3) Then
         iPL=6
      Else If (iPL.eq.4) Then
         iPL=99
      Else If (iPL.eq.5) Then
         iPL=99
      End If
      nPrint(iRout)=iPL
*
*     Reduce print level if iterating
*
      If (Reduce_Prt().and.iPL.le.5) Then
         nPrint(iRout)=4
      End If
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.6) Call CWTime(TC0,TW0)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     2 - C E N T E R   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Compute the two-center integrals over the auxiliary basis
*
      Call Drv2El_2Center_RI(Integral_RI_2,ThrAO,ipA_Diag,
     &                       nSO_Aux,MaxCntr,ipSO2C)
*
*     Post processing to generate the Q-vectors.
*
      If (LDF) Then
*
*        Local RI
*
         Call Post_2Center_LDF(ipA_Diag,ipAB,MaxCntr,Lu_AB,ipLocal_A,
     &                         nLocal_A,ipSO2C,nSO_Aux)
*
      Else
*
*        Standard RI
*
         Call Post_2Center_RI(ipA_Diag)
*
      End If
*
      Call Set_Basis_Mode('Auxiliary')
      Call Nr_Shells(nSkal_Auxiliary)
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' 2-center integrals:'
         Call CWTime(TC1,TW1)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC1-TC0,' sec.',
     &                  '      Wall time:',TW1-TW0,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     3 - C E N T E R   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      Call StatusLine(' Seward:',' Computing 3-center RI integrals')
*
*     Handle both the valence and the auxiliary basis set
*
      Call Set_Basis_Mode('WithAuxiliary')
      Call SetUp_iSD
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      Indexation = .False.
      DoGrad=.False.
      DoFock=.False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      nSkal_Valence=nSkal-nSkal_Auxiliary
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      nTMax=nSkal_Valence**2+nSkal_Auxiliary-1
      Call GetMem('TMax','Allo','Real',ipTMax,nTMax)
*
      Call Allocate_Work(ip_Tmp,nSkal**2)
      Call Shell_MxSchwz(nSkal,Work(ip_Tmp))
      TMax_all=Zero
      Do iS = 1, nSkal_Valence
         Do jS = 1, iS
            ip_Out=ip_Tmp + (jS-1)*nSkal + iS -1
            ip_In =ipTMax + (jS-1)*nSkal_Valence + iS -1
            Work(ip_In)=Work(ip_Out)
            ip_In =ipTMax + (iS-1)*nSkal_Valence + jS -1
            Work(ip_In)=Work(ip_Out)
            TMax_all=Max(TMax_all,Work(ip_Out))
         End Do
      End Do
      Do iS = 1, nSkal_Auxiliary-1
         iS_ = iS + nSkal_Valence
         jS_ = nSkal_Valence + nSkal_Auxiliary
         ip_Out = ip_Tmp + (iS_-1)*nSkal + jS_ -1
         ip_In  = ipTMax + nSkal_Valence**2 + iS -1
         Work(ip_In)=Work(ip_Out)
         TMax_all=Max(TMax_all,Work(ip_Out))
      End Do
*
      Call Free_Work(ip_Tmp)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up indexation for Gaussian pairs.
*
*     Generate some offsets and dimensions for the J12 matrix and
*     the RI vectors.
*
      Call Setup_Aux(ip_SOShl,ip_ShlSO,ip_nBasSh,nIrrep,nBas,
     &               nSkal_Valence,nSkal_Auxiliary,nSO,ip_iSSOff,
     &               Work(ipTMax),CutInt,ip_iShij,nSkal2,nBas_Aux,
     &               nChV,iTOffs)
*
      Call GetMem('iRv','Allo','Inte',ip_iRv,nSkal2)
      Call IZero(iWork(ip_iRv),nSkal2)
*                                                                      *
************************************************************************
*                                                                      *
*     Let us now decide on the memory partitioning
*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Preallocate some core for Seward!
*
      Call GetMem('MaxMem','Max','Real',iDummy,MemSew)
      MemLow=Min(MemSew/2,1024*128)
      MemSew=Max(MemSew/10,MemLow)
      Call xSetMem_Ints(MemSew)
*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     During this phase we will have three memory sections
*
*     1) the three center intergrals for a fixed {kl}
*     2) a similar block for the R-vectors
*     3) a buffer to contain subsets of the Q-vectors
*
*     Compute the max size of 1 and 2
*
      n3CMax=0
      nRvMax=0
      Call IZero(iMax_R,2*nIrrep)
      Do klS = 1, nSkal2
         kS = iWork(ip_iShij+(klS-1)*2  )
         lS = iWork(ip_iShij+(klS-1)*2+1)
         nRv = nSize_Rv(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,iOff_Rv,
     &                  nChV)
         nRvMax = Max (nRvMax,nRv)
         n3C = nSize_3C(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,iOff_3C,
     &                  nBas_Aux)
         n3CMax = Max (n3CMax,n3C)
         Do iIrrep = 0, nIrrep-1
            iMax_R(1,iIrrep)=Max(iMax_R(1,iIrrep),iOff_3C(1,iIrrep))
            iMax_R(2,iIrrep)=iMax_R(2,iIrrep)+iOff_3C(1,iIrrep)
         End Do
      End Do
*
      Call GetMem('3C','Allo','Real',ip_3C,n3CMax)
      Call GetMem('Rv','Allo','Real',ip_Rv,nRvMax)
*
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
      nQv=0
      Do iIrrep = 0, nIrrep-1
         lJ=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) lJ=lJ-1 ! remove dummy basis function
         nQv = nQv + lJ*nChV(iIrrep)
      End Do
*
*     The Q-vectors can be read in a single whole block or in chunks.
*
      If (Force_Out_of_Core) MaxMem=(8*nQv)/10
      Out_of_Core = nQv.gt.MaxMem
      nQv = Min(nQv,MaxMem)  ! note that nQv is effectively reset here
      Call GetMem('Q_vector','Allo','Real',ip_Qv,nQv)
*                                                                      *
************************************************************************
*                                                                      *
*     In case of in-core mode read Q-vectors only once!
*
      If (.Not.Out_of_Core) Then
         mQv=0
         Do iIrrep = 0, nIrrep-1
            lJ=nBas_Aux(iIrrep)
            If (iIrrep.eq.0) lJ=lJ-1 ! remove dummy basis function
*
            If (lJ.gt.0) Then
               iAddr=0
               iOff = ip_Qv+mQv
               kQv = lJ*nChV(iIrrep)
               Call dDaFile(Lu_Q(iIrrep),2,Work(iOff),kQv,iAddr)
               mQv = mQv + kQv
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Open files for the R-vectors.
*
      Do iIrrep = 0, nIrrep-1
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.ne.0) Then
            iSeed=55+iIrrep
            Lu_R(iIrrep)=IsFreeUnit(iSeed)
            Write(Name_R,'(A4,I2.2)') 'RVec',iIrrep
            Call DaName_MF_WA(Lu_R(iIrrep),Name_R)
         End If
         iAddr_R(iIrrep)=0
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*
      kCenter=0  ! dummy initialize
      lCenter=0  ! dummy initialize
      iS = nSkal ! point to dummy shell
*     Save this field for the time being!
      Call ICopy(nIrrep,iTOffs(3),3,iTtmp,1)
*
      Call Init_Tsk(id,nSkal2)
*
*
      klS=0
       iTask=0
C      Do klS = 1, nSkal2
 100     Continue
         If (.Not.Rsv_Tsk(id,klS)) Go To 200
C        Write (*,*) 'Processing shell-pair:',klS
         iTask=iTask+1
*
         iWork(ip_iRv-1+iTask) = klS
         kS = iWork(ip_iShij+(klS-1)*2  )
         lS = iWork(ip_iShij+(klS-1)*2+1)
*
*        Logic to avoid integrals with mixed muonic and electronic
*        basis.
*
         kCnttp=iSD(13,kS)
         lCnttp=iSD(13,lS)
*
         If (LDF) Then
*
*           Pick up the corresponding (K|L)^{-1} block
*
            kCenter=iSD(10,kS)
            lCenter=iSD(10,lS)
C           Write (6,*) 'kCenter, lCenter=',kCenter, lCenter
            klCenter = kCenter*(kCenter-1)/2 + lCenter
            iAdr_AB=iWork(ipAB+(klCenter-1)*2  )
            nAB    =iWork(ipAB+(klCenter-1)*2+1)
            Call dDaFile(Lu_AB,2,Work(ipLocal_A),nAB**2,iAdr_AB)
C           Call RecPrt('A^-1',' ',Work(ipLocal_A),nAB,nAB)
*
*           Now I need some lookup tables to be used below. I need to
*           go from SO index to lO index and from a given lO index
*           back to the SO index.
*
            Call IZero(ISO2LO,2*(MaxBfn+MaxBfn_Aux))
            iLO=0
            nCase=1
            If (kCenter.ne.lCenter) nCase=2
            Do iCase = 1, nCase
               If (iCase.eq.1) Then
                  jCenter=kCenter
               Else
                  jCenter=lCenter
               End If
               Do iSO_Aux = 1, nSO_Aux
                  iCenter=iWork(ipSO2C+iSO_Aux-1)
C                 Write (6,*) 'iCenter=',iCenter
                  If (iCenter.eq.jCenter) Then
                     iLO = iLO + 1
C                    Write (6,*) 'iLO,iSO_Aux=',iLO,iSO_Aux
                     iSO2LO(1,iSO_Aux)=iLO
                     iSO2LO(2,iLO)=iSO_Aux
                  End If
               End Do
            End Do
         End If
*
         Aint_kl = TMax_Valence(kS,lS)
         If (dbsc(kCnttp)%fMass.ne.dbsc(lCnttp)%fMass) Aint_kl=0.0D0
*
         nRv = nSize_Rv(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,iOff_Rv,
     &                  nChV)
         n3C = nSize_3C(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,iOff_3C,
     &                  nBas_Aux)
         Call FZero(Work(ip_3C),n3C)
         Call FZero(Work(ip_Rv),nRv)
*
         Call ICopy(nIrrep,iOff_3C,3,iTOffs(3),3)
*
*        Loop over the auxiliary basis set
*
         Do jS = nSkal_Valence+1, nSkal-1
C           Write (6,*) 'jS,kS,lS=',jS,kS,lS
            If (LDF) Then
               jCenter=iSD(10,jS)
               If (jCenter.ne.kCenter .and.
     &             jCenter.ne.lCenter ) Go To 14
C              Write (6,*) 'jCenter=',jCenter
            End If
*
            Aint=Aint_kl * TMax_Auxiliary(jS-nSkal_Valence)
*
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
            Write (6,*) 'AInt,CutInt=',AInt,CutInt
            Write (6,*)
#endif
            If (AInt.lt.CutInt) Go To 14
            Call Eval_IJKL(iS,jS,kS,lS,Work(ip_3C),n3C,Integral_WrOut)
 14         Continue
*
*           Use a time slot to save the number of tasks and shell
*           quadrupltes process by an individual node
            Call SavStat(1,One,'+')
            Call SavStat(2,One,'+')
*
         End Do    ! jS
*                                                                      *
************************************************************************
*                                                                      *
*        Multiply the 3-center integrals with the Q-vectors
*
*        Compute HQ
*
         Call Mult_3C_Qv_S(Work(ip_3C),n3C,Work(ip_Qv),nQv,Work(ip_Rv),
     &                     nRv,nChV,
     &                     iOff_3C,nIrrep,Out_of_Core,Lu_Q,'N')
*                                                                      *
************************************************************************
*                                                                      *
*        Write the R-vectors to disk. These will be retrieved and sort
*        afterwards in step 3.
*
         Do iIrrep = 0, nIrrep-1
            ip_R = ip_Rv + iOff_Rv(iIrrep)
            nRv=iOff_3C(1,iIrrep)*nChV(iIrrep)
C           Write (*,*) 'iAddr_R(iIrrep)=',iAddr_R(iIrrep)
            If (nRv.gt.0) Then
               Call dDaFile(Lu_R(iIrrep),1,Work(ip_R),nRv,
     &                      iAddr_R(iIrrep))
            End If
         End Do
*                                                                      *
************************************************************************
*                                                                      *
C      End Do    ! klS

      Go To 100
 200  Continue
      nTask=iTask
*
*     Restore iTOffs(3,*)
      Call ICopy(nIrrep,iTtmp,1,iTOffs(3),3)
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
      Call Free_Tsk(id)
*
*     Set up array to store the load balance if this might be needed in
*     a gradient calculation.
*
      Call GetMem('TmpList','Allo','Inte',ip_tmp,nSkal2)
      Call ICopy(nSkal2,[0],0,iWork(ip_tmp),1)
      Call GetMem('LBList','Allo','Inte',ip_LB,nSkal2)
      Call ICopy(nSkal2,[-1],0,iWork(ip_LB),1)
      Do iTask = 1, nTask
         klS = iWork(ip_iRv-1+iTask)
         iWork(ip_tmp-1+klS) = 1
      End Do
*
      iLB=0
      Do klS = 1, nSkal2
         If (iWork(ip_tmp-1+klS).eq.1) Then
            iWork(ip_LB+iLB)=klS
            iLB=iLB+1
         End If
      End Do
*
*
      Call Put_iArray('LBList',iWork(ip_LB),nSkal2)
*
      Call GetMem('LBList','Free','Inte',ip_LB,nSkal2)
      Call GetMem('TmpList','Free','Inte',ip_tmp,nSkal2)
*
      Call GetMem('Rv','Free','Real',ip_Rv,nRvMax)
      Call GetMem('3C','Free','Real',ip_3C,n3CMax)
      Call GetMem('Q-vector','Free','Real',ip_Qv,nQv)
      Call xRlsMem_Ints
      Call GetMem('TMax','Free','Real',ipTMax,nSkal)
      If (LDF) Then
         Call GetMem('SO2C','Free','Inte',ipSO2C,nSO_Aux)
         Call GetMem('A-blocks','Free','Inte',ipAB,
     &               MaxCntr*(MaxCntr+1)/2)
         Call GetMem('Local_A','Free','Real',ipLocal_A,2*nLocal_A)
         Call DaClos(Lu_AB)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Each node does now have an incomplete set of R-vectors!
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*
      Call Free_iWork(ip_iSSOff)
      Call Free_iWork(ip_ShlSO)
      Call Free_iWork(ip_SOShl)
      Call Free_iSD()
*
*     Let go off the Q-vectors for now!
*
      Do iIrrep = 0, nIrrep-1
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.ne.0) Call DaClos(Lu_Q(iIrrep))
      End Do
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' 3-center integrals:'
         Call CWTime(TC0,TW0)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC0-TC1,' sec.',
     &                  '      Wall time:',TW0-TW1,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     P A R T I A L   T R A N S P O S E   S E C T I O N                *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     For the interface to work fix the tables of Seward
*
      Call Set_Basis_Mode('Valence')
      Call SetUp_iSD
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do generate
*     tables for indexation.
*
      Indexation = .True.
      Call Setup_Ints(nSkal_Valence,Indexation,ThrAO,DoFock,DoGrad)
*
*     Initiate stuff for Cholesky style storage.
*
      Call IniCho_RI(nSkal_Valence,nChV,nIrrep,iTOffs,
     &               iWork(ip_iShij),nSkal2)

      Call GetMem('Addr','Allo','Inte',ip_Addr,nSkal2) ! addr for read
      Call GetMem('NuMu','Allo','INTE',iLst,2*nSkal2)
*                                                                      *
************************************************************************
*                                                                      *
*     Write out the RI vectors in Cholesky format
*
*     Here we will read one chuck from the R-vector file, while we will
*     store an as large part of the RI vectors in Cholesky format.
*
      LenVec=0
      Do iIrrep = 0, nIrrep-1
         iChoVec=0
*
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.eq.0) Go To 998
*
         iSym = iIrrep+1
*
*        NumVec: is no longer equal to the # of auxiliary functions
*
         NumVec=iTOffs(3*iIrrep+1)
         If (NumVec.eq.0) Go To 999
*
         iWork(ip_Addr) = 0
         Do i=2,nTask  ! init the addr for reading vectors
               klS = iWork(ip_iRv-1+i-1)
               kS = iWork(ip_iShij+(klS-1)*2  )
               lS = iWork(ip_iShij+(klS-1)*2+1)
               n3C = nSize_3C(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               iWork(ip_Addr+i-1) = iWork(ip_Addr+i-2) + nMuNu*NumVec
         End Do
*
         LenVec_Red = iMax_R(1,iIrrep)
         n_Rv = NumVec*LenVec_Red
         Call GetMem('Vecs','Allo','Real',ip_Rv,n_Rv)
*
*        LenVec: # of valence Gaussian products in this irrep
*
         LenVec = iMax_R(2,iIrrep)
         Call Create_Chunk(ip_iMap,ip_ChoVec,LenVec,NumVec,IncVec)
*
         Do iVec = 1, NumVec, IncVec
            NumVec_ = Min(NumVec-iVec+1,IncVec)
*                                                                      *
************************************************************************
*                                                                      *
*           Read now the R-vectors for a fixed shell-pair and
*           irrep, but for all auxiliary functions.
*
            mMuNu=0
            Do klS = 1, nSkal2
               kS = iWork(ip_iShij+(klS-1)*2  )
               lS = iWork(ip_iShij+(klS-1)*2+1)
               n3C = nSize_3C(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               m3C = nMuNu * NumVec_
*
               If (m3C.le.0) Go To 555
*
               MuNu_s=mMuNu+1
               MuNu_e=mMuNu+nMuNu
*
               iWork(iLst +(klS-1)*2  ) = MuNu_s
               iWork(iLst +(klS-1)*2+1) = MuNu_e
*
 555           Continue
               mMuNu = mMuNu + nMuNu
            End Do
*
            Do i = 1, nTask
               klS = iWork(ip_iRv-1+i)
               kS = iWork(ip_iShij+(klS-1)*2  )
               lS = iWork(ip_iShij+(klS-1)*2+1)
*
               n3C = nSize_3C(kS,lS,iWork(ip_nBasSh),nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               m3C = nMuNu * NumVec_

               If (m3C.le.0) Go To 666
*
               Call dDaFile(Lu_R(iIrrep),2,Work(ip_Rv),m3C,
     &                                       iWork(ip_Addr+i-1))

*              Copy the appropriate section into the RI vectors in
*              Cholesky format.
*
               MuNu_s = iWork(iLst +(klS-1)*2  )
               MuNu_e = iWork(iLst +(klS-1)*2+1)
               j_s=1
               j_e=NumVec_
               Call Put_Chunk(ip_ChoVec,MuNu_s,MuNu_e,
     &                        j_s,j_e,Work(ip_Rv),nMuNu,LenVec)
*
 666           Continue
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*
*           Now transfer the RI vectors to disk
*
            Call Get_Chunk(ip_ChoVec,LenVec,NumVec_,iChoVec,iSym,
     &                     ip_iMap,iVec)
*
         End Do   ! iVec = 1, NumVec, IncVec
*

         Call Destroy_Chunk(ip_ChoVec,ip_iMap)
         Call GetMem('Vecs','Free','Real',ip_Rv,n_Rv)
*
 999     Continue
*
*        Let go off the R-vectors for good!
*
         Call DaClos(Lu_R(iIrrep))
 998     Continue
         NoChoVec(iIrrep)=iChoVec
*
      End Do    ! iIrrep
      Call GetMem('NuMu','Free','INTE',iLst,2*nSkal2)
      Call GetMem('Addr','Free','Inte',ip_Addr,nSkal2)
      Call GetMem('iRv','Free','Inte',ip_iRv,nSkal2)
      Call Free_iWork(ip_nBasSh)
      Call Free_iWork(ip_iShij)
*                                                                      *
************************************************************************
*                                                                      *
      iPass = 1
      iRed = 1
      Call Cho_RI_PutInfo(iPass,iRed)
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' Block-transpose:'
         Call CWTime(TC1,TW1)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC1-TC0,' sec.',
     &                  '      Wall time:',TW1-TW0,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     D I A G O N A L   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      nDiag = 0
      Do iIrrep = 0, nIrrep-1
         nDiag = nDiag+nBas(iIrrep)
      End Do
      nDiag = nDiag*(nDiag+1)/2
      Call Allocate_Work(ipDiag,nDiag)
      Call FZero(Work(ipDiag),nDiag)
*
      Call Drv2El_RI_Diag(ThrAO,Work(ipDiag),nDiag)
*
*     Write the diagonal to disk
*
      Call Cho_IODiag(Work(ipDiag),1)
*
      Call Free_Work(ipDiag)
*
      If (iPrint.ge.6) Then
         Write (6,*) 'Diagonal vector:'
         Call CWTime(TC0,TW0)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC0-TC1,' sec.',
     &                  '      Wall time:',TW0-TW1,' sec.'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate Cholesky stuff here.
*
      irc = 0
      Call TermCho_RI(irc,NoChoVec,8)
      If (irc .ne. 0) Then
         Write(6,*) 'TermCho_RI returned ',irc
         Call SysAbendMsg('Drv2El_3Center_RI',
     &                    'Cholesky termination failed!',' ')
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
