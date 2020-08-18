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
* Copyright (C) 2007, Roland Lindh                                     *
************************************************************************

      SubRoutine Drvg1_RI(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
*     Object: superdriver for gradients for the RI/DF approximation    *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh, Dep. Chem. Phys., Lund University, Sweden  *
*             January '07                                              *
*                                                                      *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "print.fh"
#include "para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
      Real*8 Grad(nGrad), Temp(nGrad)
      Character*6 Fname
      Character*7 Fname2
      Character*8 Method
      Logical Found
      Integer nAct(0:7)
      Real*8, Allocatable:: V_k_new(:,:), U_k_new(:)
*                                                                      *
************************************************************************
*                                                                      *
      DoCholExch = .false.
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 33
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(Temp,nGrad)
      Call Allocate_Work(ipTemp,nGrad)
*                                                                      *
************************************************************************
*                                                                      *
      BufFrac=0.1D0
      Call Cho_X_Init(irc,BufFrac)
      If (irc.ne.0) Then
         Call WarningMessage(2,' Drvg1_RI: Cho_X_Init failed')
         Call Abend()
      End If
*
********************************************
*
* Decide if its MP2
*
      iMp2Prpt = 0
      Call Get_cArray('Relax Method',Method,8)
      If(Method .eq. 'MBPT2   ') Then
         Call Get_iScalar('mp2prpt',iMp2Prpt)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     In case of the Cholesky approach compute the A and Q matrices.
*
      If (Cholesky.and..Not.Do_RI) Then
*
         If (nIrrep.ne.1) Then
            Call WarningMessage(2,'Error in Drvg1_RI')
            Write (6,*) ' CD gradients with symmetry is not'
     &               //' implemented yet!'
            Call Abend()
         End If
*
         Call Cho_X_CalculateGMat(irc)
         If (iRC.ne.0) Then
            Call WarningMessage(2,'Error in Drvg1_RI')
            Write (6,*) 'Failure during G matrix construction'
            Call Abend()
         End If
*
*        Now compute the Q matrix.
*
*            Note that, as the A matrix is
*            computed in the full-pivoted (rows and columns) storage,
*            also the resulting Q matrix is full-pivoted.
*            This is necessary for the ReMap_V_k to work (see below).
*            In the RI case, only the column pivoting of Q is
*            preserved. One day we may want to unify the two cases.
*
*            (In Cholesky the Q matrix is stored as squared. In RI,
*             it is, in general, rectangular as lin. dep. may occur
*             among its columns).
*
         Call ICopy(nIrrep,NumCho,1,nBas_Aux,1)
         Call GAIGOP(nBas_Aux,nIrrep,'+')
         Call Gen_QVec(nIrrep,nBas_Aux)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Prepare handling of two-particle density.
*
      Call PrepP
*                                                                      *
************************************************************************
*                                                                      *
* Initialize the number of sets of densities and auxiliary vectors
      nAdens=1
      nAVec=1
      nKdens=1
      nJdens=1
*
      Call Qpg_iScalar('SCF mode',Found)
      If (Found) Then
         Call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
      Else
         iUHF=0
      EndIf
      nKdens=nKdens+iUHF
      nKvec=nKdens
*
      If (lPSO.and.lSA) Then
        nJdens=4
        nKdens=4
        nKVec=2
        nAdens=2
        nAvec=4
      EndIf
*
*MGD Could be more efficient memory-wise when symmetry
*--- Decompose the 2-particle active density matrix
      mAO=0
      If (lPSO) Then
         Call Get_iArray('nAsh',nAct,nIrrep)
         n_Txy=0
         Do ijsym=0,nIrrep-1
           ntmp=0
           Do jSym=0,nIrrep-1
             isym=iEOR(jSym,ijsym)
             If (iSym.gt.jSym) Then
               ntmp=ntmp+nAct(iSym)*nAct(jSym)
             Else if (iSym.eq.jSym) Then
               ntmp=ntmp+nAct(iSym)*(nAct(iSym)+1)/2
             EndIf
           EndDo
           n_Txy=n_Txy+ntmp**2
           mAO=mAO+nAct(ijsym)*nBas(ijsym)
         EndDo
         m_Txy=nAdens
         Call mma_allocate(Txy,n_Txy,nAdens,Label='Txy')
         Call mma_allocate(DMdiag,nG1,nAdens,Label='DMdiag')
         Call GetMem('Tmp','Allo','Real',ipDMtmp,nG1*(nG1+1)/2)
         Call iZero(nnP,nIrrep)
         Call Compute_txy(G1(1,1),nG1,Txy,
     &                   n_Txy,nAdens,nIrrep,DMdiag,
     &                   Work(ipDMtmp),nAct)
         Call GetMem('Tmp','Free','Real',ipDMtmp,nG1*(nG1+1)/2)
      Else
         Call mma_allocate(Txy,1,1,Label='Txy')
         Call mma_allocate(DMdiag,1,1,Label='DMdiag')
      EndIf
      n_ij2K=0
      nZ_p_k=0
      nZ_p_l=0
      nZ_p_k_New=0
      Do i=0,nIrrep-1
         iOff_ij2K(i+1) = n_ij2K
         n_ij2K = n_ij2K + nBas(i)*(nBas(i)+1)/2
         nZ_p_k = nZ_p_k + nnP(i)*nBas_Aux(i)
         nZ_p_l = nZ_p_l + nnP(i)*NumCho(i+1)
         nZ_p_k_New = nZ_p_k_New + nnP(i)*nBas(i)*(nBas(i)+1)/2
      End Do
      If (Do_RI) nZ_p_k=nZ_p_k-nnP(0)
      If (lPSO) Then
         Call mma_allocate(Z_p_k,nZ_p_k,nAvec,Label='Z_p_k')
         Z_p_k(:,:)=Zero
      Else
         Call mma_allocate(Z_p_k,1,1,Label='Z_p_k')
      EndIf
*
*     Preprocess the RI and Q vectors as follows
*
*     Allocate memory for V_k
*
      nV_k=nBas_Aux(0)
      If (Do_RI) nV_k=nV_k-1
*
      nAux_Tot=0
      Do iIrrep = 0, nIrrep-1
         nAux_Tot=nAux_Tot+nBas_Aux(iIrrep)
      End Do
*
*
      call mma_allocate(V_K,nV_k,nJdens,Label='V_k')
      V_k(:,:)=Zero
      If(iMp2prpt .eq. 2) Then
         call mma_allocate(U_K,nV_k,Label='U_k')
         U_k(:)=Zero
      Else
         call mma_allocate(U_K,1,Label='U_k')
      End If
*                    ~
*     1) Compute the V_k vector
*                 ~
*     2) Contract V_k and Q (transpose) vectors producing the V_k
*
*     Note: the above two points apply to Z_p_k as well (active space)
*
      Call GetMem('iVk','Allo','Inte',iVk,nProcs)
      Call GetMem('iZk','Allo','Inte',iZk,nProcs)
      Call IZero(iWork(iVk),nProcs)
      Call IZero(iWork(iZk),nProcs)
      iWork(iVk+myRank) = NumCho(1)*nJdens
      iWork(iZk+myRank) = nZ_p_l*nAvec
      Call GAIGOP(iWork(iVk),nProcs,'+')
      Call GAIGOP(iWork(iZk),nProcs,'+')
      iStart=1
      jStart=1
      Do j=0,nProcs-1
         itmp=iWork(iVk+j)
         iWork(iVk+j)=iStart
         iStart = iStart + itmp
         jtmp=iWork(iZk+j)
         iWork(iZk+j)=jStart
         jStart = jStart + jtmp
      End Do
*
      If(iMp2prpt .eq. 2) Then
         Call GetMem('iUk','Allo','Inte',iUk,nProcs)
         Call IZero(iWork(iUk),nProcs)
         iWork(iUk+myRank) = NumCho(1)
         Call GAIGOP(iWork(iUk),nProcs,'+')
         kStart=1
         Do j = 0,nProcs-1
            kTmp=iWork(iUk+j)
            iWork(iUk+j)=kStart
            kStart = kStart + kTmp
         End Do
      Else
         iUk = ip_iDummy
      End If
*
      Call Compute_AuxVec(iWork(iVk),iWork(iUk),iwork(iZk),
     &                    myRank+1,nProcs)
*                                                                      *
************************************************************************
*                                                                      *

      If (Cholesky.and..Not.Do_RI) Then
*
*              Map from Cholesky auxiliary basis to the full
*              1-center valence product basis.
*
         Call mma_allocate(ij2K,n_ij2K,Label='ij2K')
         ij2K(:)=0
         nV_k_New=nBas(0)*(nBas(0)+1)/2
         Call mma_allocate(V_k_new,nV_k_New,nJdens,Label="V_k_new")
         V_k_new(:,:)=Zero
*
         If (iMp2prpt .eq. 2) Then
            Call mma_allocate(U_k_new,nV_k_New,Label="U_k_new")
            U_k_new(:)=Zero
         End If

*
         Call GetMem('SO_ab','Allo','INTE',ipSO_ab,2*nAux_Tot)
         Call IZero(iWork(ipSO_ab),2*nAux_Tot)
         iOff = 0
         Do iSym = 1, nSym
            ip_List_rs=ip_InfVec+MaxVec*InfVec_N2*(iSym-1)
            Call CHO_X_GET_PARDIAG(iSym,ip_List_rs,iWork(ipSO_ab+iOff))

            If((iSym .eq. 1) .and. (iMp2prpt .eq. 2)) Then
               Call ReMap_U_k(U_k,nV_k,U_k_New,nV_k_New,
     &                        iWork(ipSO_ab))
            End If
            m_ij2K = nBas(iSym-1)*(nBas(iSym-1)+1)/2
            Do i=0,nJDens-1
              Call ReMap_V_k(iSym,V_k(1,1+i),nV_k,
     &                     V_k_new(1,1+i),nV_k_New,
     &                     iWork(ipSO_ab+iOff),ij2K(iOff_ij2K(iSym)+1),
     &                     m_ij2K)
            EndDo
            iOff = iOff + 2*nBas_Aux(iSym-1)
         End Do
*
         nV_k=nV_k_new
*
         Call Free_iWork(ipSO_ab)
         Call mma_deallocate(V_k)
         Call mma_allocate(V_k,nV_k,nJdens,Label='V_k')
         V_k(:,:)=V_k_new(:,:)
         Call mma_deallocate(V_k_new)

         If(iMp2prpt .eq. 2) Then
            Call mma_deallocate(U_k)
            Call mma_allocate(U_k,nV_k,Label='U_k')
            U_k(:)=U_k_new(:)
            Call mma_deallocate(U_k_new)
         End If

*                                                                      *
************************************************************************
*                                                                      *
*        Get the effective list of shell-pairs in case of CD
*
         Call Effective_CD_Pairs(ip_ij2,nij_Eff)
      Else
         ip_ij2 = ip_iDummy
         nij_Eff = 0
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Open C-vector-files if nSym is equal to 1
*
      If(DoCholExch) Then
         Do i=1,nKvec
           Do jSym = 1, nSym
              iSeed = 7 + jSym+(i-1)*nSym
              LuCVector(jSym,i) = IsFreeUnit(iSeed)
              If (i.eq.1) Then
                Write(Fname,'(A4,I1,I1)') 'CVEA',jSym
              ElseIf (i.eq.2) Then
                Write(Fname,'(A4,I1,I1)') 'CVEB',jSym
              EndIf
              Call DANAME_MF_WA(LuCVector(jSym,i),Fname)
            End Do
         End Do
*     Initialize timings
         Do i = 1,2
            tavec(i) = 0.0d0
            tbvec(i) = 0.0d0
         End Do
      End If
      If(imp2prpt.eq.2) Then
         Do i = 1, 2
            iSeed = 8 + nSym
            LuAVector(i) = IsFreeUnit(iSeed)
            Write(Fname2,'(A5,I1)') 'AMP2V', i
            Call DaName_MF_WA(LuAVector(i),Fname2)
            iSeed = 9 + nSym
            LuBVector(i) = IsFreeUnit(iSeed)
            Write(Fname,'(A5,I1)') 'BMP2V', i+2
            Call DaName_MF_WA(LuBVector(i),Fname)
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Compute contributions due to the "2-center" two-electron integrals
*
      Case_2C=.True.
      Call Drvg1_2center_RI(Temp,Work(ipTemp),nGrad,ip_ij2,nij_Eff)
      Call GADGOP(Work(ipTemp),nGrad,'+')
      If (iPrint.ge.15) Call PrGrad(
     &    ' RI-Two-electron contribution - 2-center term',
     &    Work(ipTemp),nGrad,lIrrep,ChDisp,iPrint)
      Call DaXpY_(nGrad,One,Temp,1,Grad,1) ! Move any 1-el contr.
      call dcopy_(nGrad,Work(ipTemp),1,Temp,1)
      Call DScal_(nGrad,-One,Temp,1)
      Case_2C=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     Compute contributions due to the "3-center" two-electron integrals
*
      Case_3C=.True.
      Call Drvg1_3center_RI(Temp,Work(ipTemp),nGrad,ip_ij2,nij_Eff)
      Call GADGOP(Work(ipTemp),nGrad,'+')
      If (iPrint.ge.15) Call PrGrad(
     &    ' RI-Two-electron contribution - 3-center term',
     &    Work(ipTemp),nGrad,lIrrep,ChDisp,iPrint)
      Call DaXpY_(nGrad,Two,Work(ipTemp),1,Temp,1)
      Case_3C=.False.
      If(Allocated(Txy))  Call mma_deallocate(Txy)
      If(Allocated(DMdiag))  Call mma_deallocate(DMdiag)
      Call GetMem('AOrb','Free','Real',ipAOrb(0,1),mAO*nADens)
*                                                                      *
************************************************************************
*                                                                      *
*
      If(DoCholExch) Then
         Do i=1,nKvec
            Do jSym = 1, nSym
               Call DaClos(luCVector(jSym,i))
            EndDo
         End Do
      End If
      If(iMp2prpt .eq. 2) Then
         Do i = 1, 2
            Call DaClos(LuAVector(i))
            Call DaClos(LuBVector(i))
         End Do
      End If


      If (Cholesky.and..Not.Do_RI) Then
         Call Free_iWork(ip_ij2)
         Call mma_deallocate(ij2K)
      End If
      Call CloseP
      Call Free_iWork(iZk)
      Call Free_iWork(iVk)
      If(iMp2prpt .eq. 2) Call Free_iWork(iUk)
      If (Allocated(Z_p_k)) Call mma_deallocate(Z_p_k)
      If (Allocated(V_k)) Call mma_deallocate(V_k)
      If (Allocated(U_k)) Call mma_deallocate(U_k)
      Call Cho_X_Final(irc)
      If (irc.ne.0) Then
         Call WarningMessage(2,' Drvg1_RI: Cho_X_Final failed')
         Call Abend()
      End If
      Call Free_Work(ipTemp)
      If (iPrint.ge.15)  Call PrGrad(
     &    ' RI-Two-electron contribution - Temp',
     &    Temp,nGrad,lIrrep,ChDisp,iPrint)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(6,TCpu2-TCpu1,TWall2-TWall1)
*
#ifdef _CD_TIMING_
      Drvg1_CPU = TCpu2-TCpu1
      Drvg1_Wall= TWall2-TWall1
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Cho_get_grad:'
      Write(6,*) 'Wall/CPU',ChoGet_Wall, ChoGet_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Mult_Rijk_Qkl:'
      Write(6,*) 'Wall/CPU',rMult_Wall, rMult_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Prepp:'
      Write(6,*) 'Wall/CPU',Prepp_Wall, Prepp_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Pget_ri2:'
      Write(6,*) 'Wall/CPU',Pget2_Wall, Pget2_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Pget_ri3:'
      Write(6,*) 'Wall/CPU',Pget3_Wall, Pget3_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Drvg1_ri:'
      Write(6,*) 'Wall/CPU',Drvg1_Wall, Drvg1_CPU
      Write(6,*) '-------------------------'
      Total_Dens_Wall = ChoGet_Wall+rMult_Wall+Prepp_Wall+Pget2_Wall +
     &                  Pget3_Wall
      Total_Dens_CPU = ChoGet_CPU+rMult_CPU+Prepp_CPU+Pget2_CPU +
     &                 Pget3_CPU
      Total_Der_Wall = Drvg1_Wall - Total_Dens_Wall
      Total_Der_CPU = Drvg1_CPU - Total_Dens_CPU
      Total_Der_Wall2 = TwoEl2_Wall + TwoEl3_Wall
      Total_Der_CPU2 = TwoEl2_CPU   + TwoEl3_CPU

      Write(6,*) 'Total Time for Density:'
      Write(6,*) 'Wall/CPU',Total_Dens_Wall, Total_Dens_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Total TIme for 2-center Derivatives:'
      Write(6,*) 'Wall/CPU',Twoel2_Wall, Twoel2_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Total TIme for 3-center Derivatives:'
      Write(6,*) 'Wall/CPU',Twoel3_Wall, Twoel3_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Total Time for Derivatives:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall2, Total_Der_CPU2
      Write(6,*) '-------------------------'
      Write(6,*) 'Derivative check:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall, Total_Der_CPU
      Write(6,*) '-------------------------'
#endif

      Return
      End
