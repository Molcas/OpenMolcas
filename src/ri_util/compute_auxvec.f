************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Compute_AuxVec(ipVk,ipZpk,myProc,nProc,ipUk)
      use pso_stuff
      use Basis_Info, only: nBas, nBas_Aux
      use Temporary_Parameters, only: force_out_of_core
      use RICD_Info, only: Do_RI, Cholesky
      use Symmetry_Info, only: nIrrep
      use Data_Structures, only: Allocate_DSBA, Map_to_DSBA
      use Data_Structures, only: Deallocate_DSBA
      use ExTerm, only: iMP2prpt
      Implicit Real*8 (a-h,o-z)
      Integer ipVk(nProc), ipZpk(nProc)
      Integer, Optional:: ipUk(nProc)
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "cholesky.fh"
#include "etwas.fh"
#include "exterm.fh"

      Type (DSBA_Type) DSQ, ChM(5)

      Integer :: ipChM(5)=[0,0,0,0,0]
      Logical DoExchange, DoCAS, Estimate, Update
      Integer nIOrb(0:7),nV_l(0:7),nV_t(0:7)
      Integer nU_l(0:7), nU_t(0:7)
      Integer ipTxy(0:7,0:7,2),ipDLT2,jp_V_k
#include "chotime.fh"
      Character*8 Method

      Real*8, Allocatable :: TmpD(:), Zv(:), Qv(:), Scr(:)
*
************************************************************************
*                                                                      *
      DoExchange=Exfac.ne.Zero
*
      nV_ls=0
      Do i=0,nIrrep-1
         nV_l(i) = NumCho(i+1) ! local # of vecs in parallel run
         nV_ls=nV_ls+nV_l(i)
         nV_t(i) = nV_l(i)
      End Do
      Call GAIGOP(nV_t,nIrrep,'+') ! total # of vecs
      If (nV_t(0).eq.0) Then
         Call WarningMessage(2,
     &          'Compute_AuxVec: no total symmetric vectors!!')
         Call Abend()
      EndIf
*
      If(iMp2prpt.eq.2) Then
          If (.NOT.Present(ipUk)) Then
            Write (6,*) 'No ipUk input present!'
            Call Abend()
          End If
          nU_ls=0
          Do i=0,nIrrep-1
             nU_l(i) = NumCho(i+1) ! local # of vecs in parallel run
             nU_ls=nU_ls+nU_l(i)
             nU_t(i) = nU_l(i)
          End Do
          Call GAIGOP(nU_t,nIrrep,'+') ! total # of vecs
          If (nU_t(0).eq.0) Then
             Call WarningMessage(2,
     &            'Compute_AuxVec: no total symmetric vectors!!')
             Call Abend()
          EndIf
      End If
*
      NChVMx=0
      nQMax=0
      Do i=0,nIrrep-1
         NChVMx= Max(NChVMx,nV_t(i))
         nQMax = Max(nQMax,nBas_Aux(i))
         nChOrb(i,1)=0
         nChOrb(i,2)=0
      End Do
      nQvMax=nQMax*NChVMx
      Call mma_allocate(Scr,nQMax)
*
      DoCAS=lPSO
*
      If (nV_ls >=1) Then ! can be = 0 in a parallel run
*
         jp_V_k = ipVk(myProc)
         jp_Z_p_k = ipZpk(myProc)
         jp_U_k = 1
         If(iMp2prpt .eq. 2) Then
            jp_U_k = ipUk(myProc)
         End If
************************************************************************
*                                                                      *
*        Get (and transform) the density matrices                      *
*                                                                      *
************************************************************************
*
         Timings=.False.
*        Timings=.True.
*
         Call Get_iArray('nIsh',nIOrb,nIrrep)

         If(iMp2prpt .ne. 2) Then
            If (DoCAS.and.lSA) Then
               nSA=5
               Call GetMem('Dens','Allo','Real',ipDMLT(1),nDens*nSA)
               Do i=2,nSA
                 ipDMLT(i)=ipDMLT(i-1)+nDens
               End Do
               call dcopy_(nDens*nSA,D0,1,Work(ipDMLT(1)),1)
*Refold some density matrices
               ij = -1
               Do iIrrep = 0, nIrrep-1
                 Do iBas = 1, nBas(iIrrep)
                   Do jBas = 1, iBas-1
                     ij = ij + 1
                     Work(ipDMLT(1)+ij)=Two*Work(ipDMLT(1)+ij)
                     Work(ipDMLT(3)+ij)=Two*Work(ipDMLT(3)+ij)
                     Work(ipDMLT(5)+ij)=Two*Work(ipDMLT(5)+ij)
                   EndDo
                   ij = ij + 1
                 EndDo
               EndDo
            Else
               Call GetMem('DMLT(1)','Allo','Real',ipDMLT(1),nDens)
               Call Get_D1AO_Var(Work(ipDMLT(1)),nDens)
            EndIf
         Else
            Call GetMem('DMLT(1)','Allo','Real',ipDMLT(1),nDens)
            Call Get_D1AO(Work(ipDMLT(1)),nDens)
         End If
*
         If (nKdens.eq.2) Then
            Call GetMem('DMLT(2)','Allo','Real',ipDMLT(2),nDens)
!           spin-density matrix
            Call Get_D1SAO_Var(Work(ipDMLT(2)),nDens)
            Call daxpy_(nDens,-One,Work(ipDMLT(1)),1,
     &                              Work(ipDMLT(2)),1)
            call dscal_(nDens,-Half,Work(ipDMLT(2)),1) ! beta DMAT
            Call daxpy_(nDens,-One,Work(ipDMLT(2)),1,
     &                             Work(ipDMLT(1)),1) ! alpha DMAT
         ElseIf (nKdens.gt.4 .or. nKdens.lt.1) Then
            Call WarningMessage(2,
     &          'Compute_AuxVec: invalid nKdens!!')
            Call Abend()
         EndIf
         ipDLT2 = 1
         If(iMp2prpt.eq.2) Then
            Call GetMem('DLT2','Allo','Real',ipDLT2,nDens)
            Call Get_D1AO_Var(Work(ipDLT2),nDens)
            Call daxpy_(nDens,-One,Work(ipDMLT(1)),1,Work(ipDLT2),1)
         Else
            ipDLT2 = ip_Dummy
         End If
************************************************************************
*                                                                      *
*       Compute Fr+In+Ac localized orbitals                            *
*       using Cholesky  decomposition for PD matrices                  *
*       using Eigenvalue decomposition for non-PD matrices (SA-CASSCF) *
*                                                                      *
************************************************************************
*         DoExchange=Exfac.ne.Zero
*
         Call Get_cArray('Relax Method',Method,8)
         If (Method.eq.'MCPDFT ' ) exfac=1.0d0
         DoExchange=Exfac.ne.Zero

         If (DoExchange .or. DoCAS) Then
            Do i=1,nKdens
              Call Allocate_DSBA(ChM(i),nBas,nBas,nSym)
              ipChM(i)=ip_of_Work(ChM(i)%A0(1))
            End Do
            If (lSA) Then
              Call mma_allocate(TmpD,nDens,Label='TmpD')
            EndIf
*                                                                      *
************************************************************************
*                                                                      *
*          PD matrices
*
            Call Allocate_DSBA(DSQ,nBas,nBas,nSym)
            Do j=1,nKvec
               If (lSA) Then
                 If (j.eq.1) Then
                    call dcopy_(nDens,Work(ipDMLT(1)),1,TmpD,1)
                 Else If (j.eq.2) Then
                    call dcopy_(nDens,Work(ipDMLT(3)),1,TmpD,1)
                 EndIf
                 Call UnFold(TmpD,nDens,DSQ%A0,SIZE(DSQ%A0),nIrrep,nBas)
               Else
                 Call UnFold(Work(ipDMLT(j)),nDens,DSQ%A0,SIZE(DSQ%A0),
     &                                      nIrrep,nBas)
               EndIf
*
               Do i=0,nIrrep-1
                  Call CD_InCore(DSQ%SB(i+1)%A2,nBas(i),
     &                            ChM(j)%SB(i+1)%A2,
     &                           nBas(i),nChOrb(i,j),1.0d-12,irc)
               End Do
               If (irc.ne.0) Then
                 Write (6,*)
     &             'Compute_AuxVec: failed to get Cholesky MOs !'
                 Call Abend()
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*          non-PD matrices
*
            If (lSA) Then
*
               Do i=3,4
*
**    Get the appropriate density matrix
*
                 If (i.eq.3) Then
                   call dcopy_(nDens,Work(ipDMLT(2)),1,TmpD,1)
                 Else If (i.eq.4) Then
                   call dcopy_(nDens,Work(ipDMLT(4)),1,TmpD,1)
                 EndIf
*
**    And eigenvalue-decompose it
*
                 iOffDSQ=0
                 Do isym=0,nIrrep-1
*
                   ChM(i)%SB(iSym+1)%A2(:,:)=Zero
                   call dcopy_(nbas(isym),[One],0,
     %                         ChM(i)%SB(iSym+1)%A2,nBas(isym)+1)
                   Call NIdiag(TmpD(1+iOffDSQ),ChM(i)%SB(iSym+1)%A2,
     &                  nBas(isym),nBas(isym),0)
*
**   First sort eigenvectors and eigenvalues
*
                   Do j=1,nBas(isym)
                     irun=iOffDSQ+j*(j+1)/2
                     Do k=j,nBas(isym)
                       jrun=iOffDSQ+k*(k+1)/2
                       If (TmpD(irun).lt.TmpD(jrun)) Then
                         tmp=TmpD(irun)
                         TmpD(irun)=TmpD(jrun)
                         TmpD(jrun)=tmp
                         Do l=1,nBas(isym)
                           tmp=ChM(i)%SB(iSym+1)%A2(l,j)
                           ChM(i)%SB(iSym+1)%A2(l,j)
     &                        = ChM(i)%SB(iSym+1)%A2(l,k)
                           ChM(i)%SB(iSym+1)%A2(l,k)=tmp
                         End Do
                       End If
                     End Do
                   End Do
*
                   Cho_thrs=1.0d-12

                   l = 0
                   Do j=1,nBas(isym)
                     If (TmpD(iOffDSQ+j*(j+1)/2).gt.Cho_thrs) Then
                       l=l+1
                       tmp=Sqrt(TmpD(iOffDSQ+j*(j+1)/2))
                       Do k=1,nBas(isym)
                         ChM(i)%SB(iSym+1)%A2(k,l) =
     &                      ChM(i)%SB(iSym+1)%A2(k,j)*tmp
                       End Do
                     EndIf
                   End Do
                   npos(isym,i-2)=l
*
                   Do j=1,nBas(isym)
                     If (-TmpD(iOffDSQ+j*(j+1)/2).gt.Cho_thrs) Then
                       l=l+1
                       irun=(l-1)*nBas(isym)
                       jrun=(j-1)*nBas(isym)
                       tmp=Sqrt(-TmpD(iOffDSQ+j*(j+1)/2))
                       Do k=1,nBas(isym)
                         ChM(i)%SB(iSym+1)%A2(k,l) =
     &                      ChM(i)%SB(iSym+1)%A2(k,j)*tmp
                       End Do
                     EndIf
                   End Do
                   nChOrb(isym,i)=l
*
                   iOffDSQ=iOffDSQ+nBas(isym)*(nBas(isym)+1)/2
                 End Do
*
               End Do
*
** Refold the other DM
*
               ij = -1
               Do iIrrep = 0, nIrrep-1
                 Do iBas = 1, nBas(iIrrep)
                   Do jBas = 1, iBas-1
                     ij = ij + 1
                     Work(ipDMLT(2)+ij)=Two*Work(ipDMLT(2)+ij)
                     Work(ipDMLT(4)+ij)=Two*Work(ipDMLT(4)+ij)
                   EndDo
                   ij = ij + 1
                 EndDo
               EndDo
            EndIf
            Call Deallocate_DSBA(DSQ)
            If (lSA) Call mma_deallocate(TmpD)
         EndIf
************************************************************************
*                                                                      *
*        First contract the RI vectors with the density matrix         *
*                                                                      *
************************************************************************
* --- Pointers to the Cholesky vectors of P2
         mAO=0
         iOff=0
         Do kIrrep=0,nIrrep-1 ! compound symmetry
            iOff2=0
            Do jIrrep=0,nIrrep-1
               iIrrep=iEOR(jIrrep,kIrrep)
               If (iIrrep.lt.jIrrep) Then
                  nnAorb=nASh(iIrrep)*nAsh(jIrrep)
               ElseIf (iIrrep.eq.jIrrep) Then
                  nnAorb=nAsh(iIrrep)*(nAsh(iIrrep)+1)/2
               Else
                  Go To 100
               EndIf
               ipTxy(iIrrep,jIrrep,1) = 1 + iOff2+iOff
               ipTxy(jIrrep,iIrrep,1) = ipTxy(iIrrep,jIrrep,1)
               If (lSA) Then
                 ipTxy(iIrrep,jIrrep,2) = ipTxy(iIrrep,jIrrep,1)+n_Txy
                 ipTxy(jIrrep,iIrrep,2) = ipTxy(iIrrep,jIrrep,2)
               EndIf
               iOff2=iOff2+nnAorb
100            Continue
            End Do
            iOff=iOff+iOff2*nnP(kIrrep)
            mAO=mAO+nBas(kIrrep)*nASh(kIrrep)
         End Do

         Allocate(AOrb(nADens))
         Do iADens = 1, nADens
            Call Allocate_DSBA(AOrb(iADens),nAsh,nBas,nIrrep)
            Call Map_to_DSBA(AOrb(iADens),ipAOrb(:,iADens))
         End Do

*
* --- Reordering of the active MOs :  C(a,v) ---> C(v,a)
*
         iCount=0
         Do iIrrep=0,nIrrep-1

            jCount = iCount + nBas(iIrrep)*nIOrb(iIrrep)

            iCount = iCount + nBas(iIrrep)**2
            If (nBas(iIrrep)*nASh(iIrrep)==0) Cycle

         Do iADens=1, nADens

            Do i=1,nASh(iIrrep)
               kOff1 = 1 + jCount + nBas(iIrrep)*(i-1)
               AOrb(iADens)%SB(iIrrep+1)%A2(i,1:nBas(iIrrep)) =
     &            CMO(kOff1:kOff1-1+nBas(iIrrep),iADens)
            End Do

         End Do ! iADens

         End Do ! iIrrep
*
         If (nKdens.eq.2) Then ! for Coulomb term
            Call daxpy_(nDens,One,Work(ipDMLT(2)),1,
     &                            Work(ipDMLT(1)),1)
         EndIf
*
*  Add the density of the environment in a OFE calculation (Coulomb)
*
         Call OFembed_dmat(Work(ipDMlt(1)),nDens)
*
*        nScreen=10 ! Some default values for the screening parameters
*        dmpK=One
         Estimate=.False.
         Update=.True.
         Call Cho_Get_Grad(irc,nKdens,ipDMlt,ipDLT2,ipChM,
     &                     Txy,n_Txy*nAdens,ipTxy,
     &                     DoExchange,lSA,nChOrb,AOrb,nAsh,
     &                     DoCAS,Estimate,Update,
     &                     V_k(jp_V_k,1), nV_k,
     &                     U_k(jp_U_k),
     &                     Z_p_k(jp_Z_p_k,1), nZ_p_k,
     &                     nnP,npos)
*
         If (irc.ne.0) Then
            Call WarningMessage(2,
     &                  'Compute_AuxVec: failed in Cho_Get_Grad !!')
            Call Abend()
         End If

         If (DoCAS .or. DoExchange) Then
            Do i=1,nKdens
              Call deallocate_DSBA(ChM(i))
            End Do
         EndIf
         If(iMp2prpt.eq.2) Then
            Call Free_Work(ipDLT2)
         End If
*
      End If ! no vectors on this node
*
*     For parallel run: reordering of the V_k(tilde) vector from
*     the "node storage" to the Q-vector storage
      If (nProc.gt.1)  Then
         Do i = 1, SIZE(V_K,2)
            Call Reord_Vk(ipVk,nProc,myProc,nV_l,nV_t,[1],1,V_k(:,i))
         End Do
      End If
************************************************************************
*                                                                      *
*     Second step: contract with the Q-vectors to produce V_k          *
*            ~  T                                                      *
*     V = V Q                                                          *
*                                                                      *
************************************************************************
*
      If (Cholesky.and..Not.Do_RI) Then ! to cope with the calls below
         nBas_Aux(0)=nBas_Aux(0)+1
      End If
*
      Call mma_maxDBLE(MemMax)
*
      If (Force_out_of_Core) MemMax=4*(nQvMax)/10
      nQv = Min(MemMax,nQvMax)
      Call mma_allocate(Qv,nQv,Label='Qv')
*
**    Coulomb
*
      Do i=0,nJdens-1
         Call Mult_Vk_Qv_s(V_k(ipVk(1),1+i),nV_t(0),Qv,nQv,
     &                     Scr,nQMax,nBas_Aux,nV_t(0),nIrrep,'T')
         call dcopy_(nV_k,Scr,1,V_k(ipVk(1),1+i),1)
      End Do
*
**    MP2
*
      If(iMp2prpt.eq.2) Then
         Call Mult_Vk_Qv_s(U_k(ipUk(1)),nU_t(0),Qv,nQv,
     &                     Scr,nQMax,nBas_Aux,nU_t(0),nIrrep,'T')
         call dcopy_(nV_k,Scr,1,U_k(ipUk(1)),1)
      End If
*
**    Active term
*
      If (DoCAS) Then ! reorder Zp_k

         Call mma_allocate(Zv,nZ_p_k,Label='Zv')
*
         Do iAvec=1,nAvec
           If (nProc.gt.1) Call Reord_Vk(ipZpk(1),nProc,myProc,
     &                    nV_l,nV_t,nnP,nIrrep,Z_p_k(:,iAVec))
*
           Call Mult_Zp_Qv_s(Z_p_k(ipZpk(1),iAvec),nZ_p_k,Qv,nQv,
     &                       Zv,nZ_p_k,nV_t,nnP,nBas_Aux,nIrrep,'T')
*
          call dcopy_(nZ_p_k,Zv,1,Z_p_k(ipZpk(1),iAvec),1)
         End Do
         Call mma_deallocate(Zv)
      EndIf
*
      Call mma_deallocate(Qv)
      Call mma_deallocate(Scr)
*
**    Exchange
*
      If (DoExchange) Then
         DoCholExch = .true.
         Do iSO=1,nKvec
            Call Mult_RijK_QKL(iSO,nBas_aux,nIrrep)
         End Do
         If(iMp2prpt.eq.2) Then
            Call Mult_with_Q_MP2(nBas_aux,nBas,nIrrep)
         End If
      End If
      If (Cholesky.and..Not.Do_RI) nBas_Aux(0)=nBas_Aux(0)-1
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
      Subroutine OFembed_dmat(Dens,nDens)

      use OFembed, only: Do_OFemb
      Implicit Real*8 (a-h,o-z)
      Real*8 Dens(nDens)
#include "stdalloc.fh"
      Character*16 NamRfil
      Real*8, Allocatable :: D_Var(:)

      If (.not.Do_OFemb) Return
*
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name

      Call mma_allocate(D_var,nDens,Label='D_var')
      Call get_dArray('D1aoVar',D_var,nDens)
      Call daxpy_(nDens,One,D_var,1,Dens,1)
      Call mma_deallocate(D_Var)
*
      Call NameRun(NamRfil)   ! switch back to old RUNFILE name
*
      Return
      End
