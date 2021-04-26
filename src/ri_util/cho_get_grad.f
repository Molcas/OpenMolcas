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
* Copyright (C) 2007, Francesco Aquilante                              *
*               2011, Thomas Bondo Pedersen                            *
************************************************************************
      SUBROUTINE CHO_GET_GRAD(irc,nDen,
     &                        DLT,DLT2,MSQ,
     &                        Txy,nTxy,ipTxy,DoExchange,lSA,
     &                        nChOrb_,AOrb,nAorb,DoCAS,
     &                        Estimate,Update,
     &                        V_k,nV_k,
     &                        U_k,
     &                        Z_p_k,nZ_p_k,
     &                        nnP,npos)

************************************************************************
*  Author : F. Aquilante (visiting F. Illas group in Barcelona, Spain, *
*                                                    March-April 2007) *
*                                                                      *
*  Purpose:                                                            *
*         Computation of the relevant quantities for RI                *
*         (and Cholesky) gradient code                                 *
*                                                                      *
*         Coulomb term :  V_k = sum_gd  D_gd L_gd_k                    *
*                                                                      *
*         MP2 Coulomb term : U_k = sum_gd D(MP2)_gd L_gd_k             *
*                                                                      *
*         Active term :   Z_p_k = sum_xy T(xy,p) L_xy_k                *
*                                                                      *
*         Inact. Exchange term: the quantity returned on disk is       *
*                                                                      *
*                     L_ij_k = sum_gd L_gd_k C_gi C_dj                 *
*                                                                      *
*                                                                      *
*  Input:                                                              *
*                                                                      *
*         nDen : is equal to 2 iff Spin Unrestricted                   *
*                4 for SA-CASSCF, otherwise nDen=1                     *
*                                                                      *
*         DLT : the LT-packed and symm. blocked one-body Dmat.         *
*                 For spin unrestricted, Dmat = Dalpha + Dbeta         *
*                                                                      *
*         DLT2: pointer to the LT-packed and symm. blocked             *
*                 one body MP2 Dmat.                                   *
*                                                                      *
*         MSQ : Cholesky MOs stored as C(a,i) symm. blocked            *
*                 with blocks of dim. (nBas,nBas). These are           *
*                 obtained from CD of the 1-particle DMAT.             *
*                 (Two pointers iff alpha and beta spinorbitals)       *
*                                                                      *
*         ipTxy : array (8,8) of pointers to the symm. blocks          *
*                 of the Cholesky decomposed MO-basis (symmetrized)    *
*                 2-body density matrix                                *
*                 T(xy,p) : is stored by compound symmetry JSYM        *
*                            the indices {xy} are stored as PACKED     *
*                            (sym x.le.sym y)                          *
*                                                                      *
*         DoExchange : logical to activate exchange grad. components   *
*                                                                      *
*         nChOrb_ : array of nr. of Cholesky orbitals in each irrep    *
*                                                                      *
*         nAorb : array with # of Active orbitals in each irrep        *
*                 (The same orbital basis                              *
*                 in which the 2-body Dmat is expressed)               *
*                                                                      *
*         DoCAS : logical to activate CASSCF grad. components          *
*                                                                      *
*         nScreen : See e.g. LK-screening docum. in SCF                *
*                   or CASSCF read-input routines. Default = 10        *
*                                                                      *
*         dmpK : damping for the LK-screening threshold. Def: 1.0d0    *
*                                                                      *
*         Estimate : logical for LK-screening. Default: .false.        *
*                                                                      *
*         Update : logical for LK-screening. Default: .true.           *
*                                                                      *
*         nnP : array of # of Cholesky vectors for the dec 2-body      *
*               density matrix in each compound symmetry               *
*                                                                      *
*                                                                      *
*  Output:                                                             *
*         irc : return code                                            *
*                                                                      *
*         V_k : array Real*8 for the Coulomb interm. Size=NumCho(1)    *
*                                                                      *
*         U_k : array Real*8 for the mp2 Coulomb interm. Size=NumCho(1)*
*                                                                      *
*         Z_p_k : array Real*8 for the active grad. components.        *
*                  Must be zeroed by the calling routine. Stored       *
*                  according to jSym and blocked after symm. blocks    *
*                  of the active orbitals (square storage).            *
*                                                                      *
*  Modifications:                                                      *
*    August 24, 2011, Thomas Bondo Pedersen:                           *
*       Allow zero vectors on a node.                                  *
*                                                                      *
************************************************************************
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: nnBstRSh, InfVec, IndRed
      use Data_Structures, only: DSBA_Type, SBA_Type
      use Data_Structures, only: Allocate_SBA, Deallocate_SBA
      use Data_Structures, only: Allocate_DSBA
      use Data_Structures, only: NDSBA_Type, Allocate_NDSBA,
     &                           Deallocate_NDSBA
      use Data_Structures, only: Allocate_L_Full, Deallocate_L_Full,
     &                           L_Full_Type
      use Data_Structures, only: Allocate_Lab, Deallocate_Lab,
     &                           Lab_Type
      use ExTerm, only: VJ, iMP2prpt, CMOi
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)

      Type (NDSBA_Type) DiaH
      Type (DSBA_Type) DLT(5), DLT2, MSQ(nDen), AOrb(*)
      Type (SBA_Type) Laq(1), Lxy
      Type (L_Full_Type) L_Full
      Type (Lab_Type) Lab

      Logical   DoExchange,DoCAS,lSA
      Logical   DoScreen,Estimate,Update,BatchWarn
      Integer   nDen,nChOrb_(8,5),nAorb(8),nnP(8),nIt(5)
      Integer   ipTxy(8,8,2)
      Integer   kOff(8,5), LuRVec(8,3)
      Integer :: ipDLT(5)=[0,0,0,0,0],ipDLT2=0
      Integer   npos(8,3)
      Integer   iSTSQ(8), nnA(8,8), nInd
      Real*8    tread(2),tcoul(2),tmotr(2),tscrn(2),tcasg(2),tmotr2(2)

      Real*8    Txy(nTxy),V_k(nV_k,*),Z_p_k(nZ_p_k,*), U_k(*)

      Character*6  Fname
      Character*50 CFmt
      Character(LEN=12), Parameter :: SECNAM = 'CHO_GET_GRAD'
#include "chotime.fh"
#include "real.fh"

      Logical, Parameter :: DoRead = .false.
      Real*8, Parameter :: xone = -One
#include "itmax.fh"
#include "Molcas.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"
#include "exterm.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "print.fh"
#include "bdshell.fh"
      Logical add
      Character(LEN=6) mode

      Real*8, Allocatable:: Lrs(:,:), Drs(:,:), Diag(:), AbsC(:),
     &                      SvShp(:,:), MLk(:), Ylk(:,:), Drs2(:,:)
      Real*8, Allocatable, Target:: Yik(:)
      Real*8, Pointer:: pYik(:,:)=>Null()
      Integer, Allocatable:: kOffSh(:,:), iShp_rs(:),
     &                       Indx(:,:), Indik(:,:)
      Real*8, Allocatable, Target:: Aux(:)
      Real*8, Pointer:: Lik(:,:), Rik(:)
#if defined (_MOLCAS_MPP_)
      Real*8, Allocatable:: DiagJ(:)
#endif

      Type V2
        Real*8, Pointer :: A2(:,:)=>Null()
      End Type V2

      Type Special
        Real*8, Allocatable:: A0(:)
        Type (V2) :: Den(5)
      End Type Special

      Type (Special), Target:: SumClk
*                                                                      *
************************************************************************
*                                                                      *
      Interface

        Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                         iSwap,IREDC,nDen,kDen,MOs,ChoT,
     &                         DoRead)
        use Data_Structures, only: DSBA_Type, SBA_Type
        Integer irc, lRedVec
        Real*8 RedVec(lRedVec)
        Integer IVEC1,NUMV,ISYM,iSwap,IREDC
        Integer   nDen,kDen

        Type (DSBA_Type) MOs(nDen)
        Type (SBA_Type) Chot(nDen)

        Logical   DoRead
        End Subroutine Cho_X_getVtra

        SUBROUTINE CHO_GetShFull(LabJ,lLabJ,JNUM,JSYM,IREDC,ChoV,
     &                          SvShp,mmShl,iShp_rs,mmShl_tot)
        use Data_Structures, only: L_Full_Type
        Integer lLabJ, JNUM, JSYM, IREDC
        Integer mmShl, mmShl_tot
        Real*8  LabJ(lLabJ)
        Type (L_Full_Type) ChoV
        Real*8  SvShp(mmShl , 2)
        Integer iShp_rs( mmShl_tot )
        End SUBROUTINE CHO_GetShFull

      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1

      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*                                                                      *
************************************************************************
*                                                                      *
*     General Initialization                                           *
*                                                                      *
************************************************************************

      Do i = 1, 5
         If (DLT(i)%Active) ipDLT(i)=ip_of_Work(DLT(i)%A0(1))
      End Do
      If (DLT2%Active) ipDLT2=ip_of_Work(DLT2%A0(1))

      iRout = 9
      iPrint = nPrint(iRout)

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read vectors
      tcoul(:) = zero  !time for computing V_k
      tcasg(:) = zero  !time for computing Z_p_k
      tmotr(:) = zero  !time for the MO transf of vectors
      tmotr2(:)= zero  !time for the 2nd MO transf of vectors
      tscrn(:) = zero  !time for screening overhead

      IREDC = -1  ! unknown reduced set in core

      BatchWarn = .True.
      nInd = 0

      Call set_nnA(nSym,nAorb,nnA)
*
**    Various offsets
*
      MaxB=nBas(1)
      ISTSQ(1)=0
      DO ISYM=2,NSYM
        MaxB=Max(MaxB,nBas(iSym))
        NBQ=NBAS(ISYM-1)**2
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NBQ ! Diagonal integrals in full
      END DO
*
**
*
      nI2t=0
      nItmx=0
      nIt(:) = 0
      Do jDen=nDen,1,-1
         kOff(1,jDen)=0
         nIt(jDen)=nChOrb_(1,jDen)
         Do i=2,nSym
            kOff(i,jDen)=nIt(jDen)
            nIt(jDen)=nIt(jDen)+nChOrb_(i,jDen)
         End Do
         nI2t=nI2t+nIt(jDen)
         nItmx=Max(nItmx,nIt(jDen))
      End Do
*
**   Initialize pointers to avoid compiler warnings
*
      thrv=0.0d0
      xtau=0.0d0
*
**    Construct iBDsh for later use
*
      Do iSyma=1,nSym
         LKsh=0
         Do iaSh=1,nShell
            iSSa=nShell*(iSyma-1)+iaSh
            iBDsh(iSSa) = LKsh
            LKsh = LKsh + nBasSh(iSyma,iaSh)
         End Do
      End Do

!     iShp_rs
      Call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

************************************************************************
*                                                                      *
*     Initialize a few things for ij-screening //Jonas B               *
*                                                                      *
************************************************************************
      If(DoExchange) Then
*
** Define the screening thresholds
*

         Call Get_dScalar('Cholesky Threshold',ThrCom)

         tau = (ThrCom/DBLE(Max(1,nItmx)))*dmpK

         MaxRedT=MaxRed
         Call GAIGOP_SCAL(MaxRedT,'+')

         If (Estimate) tau=tau/DBLE(MaxRedT)
         xtau = sqrt(tau)

         NumVT=NumChT
         Call GAIGOP_SCAL(NumVT,'+')
!        Vector MO transformation screening thresholds
         thrv = ( sqrt(ThrCom/DBLE(Max(1,nItmx)*NumVT)) )*dmpK

#if defined (_MOLCAS_MPP_)
         If (Is_Real_Par() .and. Update) Then
            NNBSTMX=0
            Do i=1,nSym
               NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
            End Do
            Call mma_allocate(DiagJ,NNBSTMX,Label='DiagJ')
            DiagJ(:)=Zero
         EndIf
#endif

*
** Read the diagonal integrals (stored as 1st red set)
*
         Call mma_allocate(DIAG,NNBSTRT(1),Label='Diag')
         If (Update) CALL CHO_IODIAG(DIAG,2) ! 2 means "read"

*
** Allocate memory
*
!        sqrt(D(a,b)) stored in full (squared) dim
         Call Allocate_NDSBA(DiaH,nBas,nBas,nSym)
         DiaH%A0(:)=Zero

         Call mma_allocate(AbsC,MaxB,Label='AbsC')

         Call mma_allocate(Ylk,MaxB,nItmx,Label='Ylk')

         Call mma_allocate(Yik,nItmx**2,Label='Yik') ! Yi[k] vectors

*used to be nShell*something
!        ML[k] lists of largest elements in significant shells
         Call mma_allocate(MLk,nShell,Label='MLk')

!        list of S:= sum_l abs(C(l)[k])
         Call mma_allocate(SumClk%A0,nShell*nI2t,Label='SumClk%A0')
         iE = 0
         Do i=1,nDen
           iS = iE + 1
           iE = iE + nShell*nIt(i)
           SumClk%Den(i)%A2(1:nShell,1:nIt(i))
     &          => SumClk%A0(iS:iE)
         End Do

*
** Indx and Indik must be stored for each density, symmetry, etc.
** in case of a batched procedure
*
         Do jDen=1,nKvec
           Do kSym=1,nSym
             nInd = nInd + nChOrb_(kSym,jDen)
           End Do
         End Do

!        Index array
         Call mma_allocate(Indx,[0,nShell],[1,nInd],Label='Indx')

         !Yi[k] Index array
         Call mma_allocate(Indik,(nItmx+1)*nItmx+1,nInd,Label='Indik')

!        kOffSh
         Call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

!        shell-pair Frobenius norm of the vectors
         Call mma_allocate(SvShp,nnShl,2,Label='SvShp')

*
** Jonas - June 2010:
** allocate memory for rearranged CMO-matrix
*
         Do i=1,nDen
           Call Allocate_DSBA(CMOi(i),nChOrb_(:,i),nBas,nSym)
         End Do

         nQoT = 0
*
** Compute Shell Offsets ( MOs and transformed vectors)
*
         Do iSyma=1,nSym
            LKsh=0
            Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

               kOffSh(iaSh,iSyma) = LKsh

               LKsh = LKsh + nBasSh(iSyma,iaSh)
            End Do
         End Do

*
** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
*
         Do jDen=1,nDen
            Do kSym=1,nSym

               Do jK=1,nChOrb_(kSym,jDen)
                  jK_a = jK + kOff(kSym,jDen)
*
                  Do iaSh=1,nShell

                     iS = kOffSh(iaSh,kSym) + 1
                     iE = kOffSh(iaSh,kSym) + nBasSh(kSym,iaSh)

                     SKSh=Zero
                     Do ik=iS,iE
                        SKsh = SKsh
     &                       + MSQ(jDen)%SB(kSym)%A2(ik,jK)**2
                     End Do

                     SumClk%Den(jDen)%A2(iaSh,jK_a) = SkSh

                  End Do
               End Do
            End Do
         End Do
*
** Reorder CMO-matrix, Needed to construct B-matrix for exchange
** Jonas - June 2010
*
         Do jDen = 1, nKdens
            Do kSym = 1, nSym
*
**If the orbitals come from eigenvalue decomposition, change sign
*
               If (lSA.and.(jDen.ge.3)) Then
                 npos2=npos(ksym,jDen-2)
                 Do jK = 1, nPos2
                    Do jGam = 1, nBas(kSym)
                       CMOi(jDen)%SB(kSym)%A2(jK,jGam) =
     &                               MSQ(jDen)%SB(kSym)%A2(jGam,jK)
                    End Do
                 End Do
                 Do jK = npos2+1,nChOrb_(kSym,jDen)
                    Do jGam = 1, nBas(kSym)
                       CMOi(jDen)%SB(kSym)%A2(jK,jGam) =
     &                              - MSQ(jDen)%SB(kSym)%A2(jGam,jK)
                    End Do
                 End Do
               Else
*
                 Do jK = 1, nChOrb_(kSym,jDen)
                    Do jGam = 1, nBas(kSym)
                       CMOi(jDen)%SB(kSym)%A2(jK,jGam) =
     &                               MSQ(jDen)%SB(kSym)%A2(jGam,jK)
                    End Do
                 End Do
               EndIf
            End Do
         End Do
      End If
*
** Mapping shell pairs from the full to the reduced set
*
      Call Mk_iShp_rs(iShp_rs,nShell)

************************************************************************
*                                                                      *
*     BIG LOOP OVER VECTORS SYMMETRY                                   *
*                                                                      *
************************************************************************
*                                                                      *
      DO jSym=1,nSym
*                                                                      *
************************************************************************
*                                                                      *
         NumCV=NumCho(jSym)
         Call GAIGOP_SCAL(NumCV,'max')
         If (NumCV .lt. 1) Cycle
*
** offsets for active term
*
         iOffZp=0
         Do j=1,jSym-1
            iOffZp = iOffZp + nnP(j)*NumCho(j)
         End Do
*
** Open some files to store exchange auxiliary vectors
*
         If (DoExchange) Then
            iSeed=7
            Do i=1,nSym
               k=muld2h(jSym,i)
               LuRVec(i,1) = IsFreeUnit(iSeed)
               Write(Fname,'(A4,I1,I1)') 'CHTA',i,k
               Call DANAME_MF_WA(LuRVec(i,1),Fname)
               iSeed=iSeed+1
               If (nKvec.ge.2) Then
                  LuRVec(i,2) = IsFreeUnit(iSeed)
                  Write(Fname,'(A4,I1,I1)') 'CHTB',i,k
                  Call DANAME_MF_WA(LuRVec(i,2),Fname)
                  iSeed=iSeed+1
               EndIf
            Enddo
         EndIf
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*          M E M O R Y   M A N A G E M E N T   S E C T I O N           *
*                                                                      *
************************************************************************
************************************************************************
*
*        For one Cholesky vector, JNUM=1, compute the amount of memory
*        needed for the various vectors.

         JNUM=1

         ! L_Full
         Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,
     &                        Memory=nL_Full)
         ! Lab
         mDen=1
         Call Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen,
     &                        Memory=nLab)
         If (DoCas) Then
            iSwap = 0  ! Lvb,J are returned
            Call Allocate_SBA(Laq(1),nAorb,nBas,JNUM,JSYM,nSym,
     &                        iSwap,Memory=nLaq)
            nLxy=0
            Do iMO1=1,nAdens
               iSwap_lxy=5
               If (iMO1==2) iSwap_lxy=6
               Call Allocate_SBA(Lxy,nAorb,nAorb,JNUM,JSYM,nSym,
     &                           iSwap_lxy,Memory=nLxy0)
               nLxy = Max( nLxy, nLxy0)
           End Do
         Else
            nLaq=0
            nLxy=0
         End If
*
** compute memory needed to store at least 1 vector of JSYM
** and do all the subsequent calculations
*
         nLik=0
         nRik=0
         do l=1,nSym
            k=Muld2h(l,JSYM)
            Do jDen=1,nDen
                nRik=Max(nRik,nChOrb_(l,jDen)*nChOrb_(k,jDen))
                If (nChOrb_(k,jDen).gt.0) Then
                   nLik=Max(nLik,nChOrb_(l,jDen))
                EndIf
            End Do
         end do

         ! re-use memory for the active vec
         LFMAX = Max(   nLaq + nLxy,  nL_Full + nRik + nLik + nLab )
*                                                                      *
************************************************************************
************************************************************************
*
*
**
*
         iLoc = 3 ! use scratch location in reduced index arrays

         If (NumCho(jSym).lt.1) Then
            JRED1 = 1
            JRED2 = 1
         Else
            JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
            JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last
*                                               !vec
         End If
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif

c ---    entire red sets range for parallel run
         Call GAIGOP_SCAL(JRED1,'min')
         Call GAIGOP_SCAL(JRED2,'max')
*
**       MGD does it need to be so?
*
         DoScreen=.True.
         kscreen=1
*                                                                      *
************************************************************************
*                                                                      *
         Do JRED=JRED1,JRED2
*                                                                      *
************************************************************************
*                                                                      *

            If (NumCho(jSym).lt.1) Then
               iVrs=0
               nVrs=0
            Else
               CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)
            End If

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//
     &          ': Cho_X_nVecRS returned nVrs<0. STOP!'
               call Abend()
            endif

c           !set index arrays at iLoc
            Call Cho_X_SetRed(irc,iLoc,JRED)
            if(irc.ne.0)then
              Write(6,*) SECNAM,': cho_X_setred non-zero return code.',
     &                   ' rc= ',irc
              call Abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then
               Call mma_allocate(Drs,nRS,nJdens,Label='Drs')
               Drs(:,:)=Zero
               If(iMp2prpt.eq.2) Then
                  Call mma_allocate(Drs2,nRS,1,Label='Drs2')
               End If
            End If

            Call mma_maxDBLE(LWORK)

            nVec = Min(LWORK/(nRS+LFMAX),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) ' LWORK= ',LWORK
               WRITE(6,*) ' min. mem. need= ',nRS+LFMAX
               WRITE(6,*) ' jsym= ',jsym
               WRITE(6,*) ' nRS = ',nRS
               WRITE(6,*) ' LFMAX = ',LFMAX
               WRITE(6,*)
               WRITE(6,*) ' nL_Full = ',nL_Full
               WRITE(6,*) ' nRik = ',nRik
               WRITE(6,*) ' nLik = ',nLik
               WRITE(6,*) ' nLab = ',nLab
               WRITE(6,*)
               WRITE(6,*) ' nLaq = ',nLaq
               WRITE(6,*) ' nLxy = ',nLxy
               irc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

*                                                                      *
************************************************************************
*                                                                      *
            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
C --- Transform the densities to reduced set storage
               mode = 'toreds'
               add  = .false.
               nMat=1
               Do jDen=1,nJdens
                  Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                              [ipDLT(jDen)],Drs(:,jDen),
     &                              mode,add)
               End Do
               If(iMp2prpt .eq. 2) Then
                  Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                              [ipDLT2],Drs2(:,1),mode,add)
               End If
            EndIf
*
**  BATCH over the vectors
*

            nBatch = (nVrs-1)/nVec + 1

            If (BatchWarn .and. nBatch.gt.1) Then
               If (iPrint.ge.6) Then
                  Write(6,'(20A3)')('---',I=1,20)
                  Write(6,*)' Batch procedure used.'//
     &                      ' Increase memory if possible!'
                  Write(6,'(20A3)')('---',I=1,20)
                  Write(6,*)
                  Call XFlush(6)
               End If
               BatchWarn = .False.
            EndIf

*                                                                      *
************************************************************************
*                                                                      *
            DO iBatch=1,nBatch
*                                                                      *
************************************************************************
*                                                                      *
               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               Else
                  JNUM = nVec
               Endif

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  irc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

************************************************************************
************************************************************************
**                                                                    **
**                                                                    **
**             Coulomb term                                           **
**             V{#J} = sum_ab  L(ab,{#J}) * D(ab)                     **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
               If(JSYM.eq.1)Then

                 CALL CWTIME(TCC1,TWC1)

*
**  Inactive Coulomb term
*
                 Do jden=1,nJdens
                    CALL DGEMV_('T',nRS,JNUM,
     &                         One,Lrs,nRS,
     &                             Drs(1,jden),1,
     &                         zero,V_k(jVec,jDen),1)
                 End Do
*
**  MP2 Coulomb term
*
                 If(iMp2prpt .eq. 2) Then
                    CALL DGEMV_('T',nRS,JNUM,
     &                         One,Lrs,nRS,
     &                             Drs2(:,1),1,
     &                         zero,U_k(jVec),1)
                 End If
*
                 CALL CWTIME(TCC2,TWC2)
                 tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                 tcoul(2) = tcoul(2) + (TWC2 - TWC1)
               EndIf
************************************************************************
************************************************************************
**                                                                    **
**             E X C H A N G E    T E R M                             **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
*
               If (DoExchange) Then

                  CALL CWTIME(TCS1,TWS1)
************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      Select only important  ij pairs                                 *
*      For this, one computes the quantity                             *
*         Yik = sum_mu_nu (mu nu | mu nu)^1/2 X_mu_i X_nu_k            *
*      with (mu nu | mu nu) = sum_J (L_mu_nu,J)^2                      *
*                                                                      *
*                                                                      *
*      a) Estimate the diagonals :   D(mu,nu) = sum_J (L_mu_nu,J)^2    *
*                                                                      *
************************************************************************
                  If (Estimate) Then

                     Call Fzero(Diag(1+iiBstR(jSym,1)),
     &                          NNBSTR(jSym,1))

                     Do krs=1,nRS

                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set

                        Do jvc=1,JNUM

                           Diag(jrs) = Diag(jrs) + Lrs(krs,jvc)**2

                        End Do

                     End Do

                  EndIf

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
                  Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,
     &                                 nSym)
                  Call mma_allocate(Aux,(nRik+nLik)*nVec,Label='Aux')
                  Call Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,
     &                              mDen)

                  CALL CWTIME(TCX1,TWX1)

*
** Reorder vectors to Full-dimensions
**
** Vectors are returned in the storage LaJ,b with the restriction:
**    Sym(a).ge.Sym(b)
** and blocked in shell pairs
*
                  CALL CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,
     &                               SvShp,nnShl,iShp_rs,nnShl_tot)

                  CALL CWTIME(TCX2,TWX2)
                  tmotr(1) = tmotr(1) + (TCX2 - TCX1)
                  tmotr(2) = tmotr(2) + (TWX2 - TWX1)

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      b) DH(mu,nu)=sqrt(D(mu,nu))                                     *
*         Only the symmetry blocks with compound symmetry JSYM         *
*         are computed                                                 *
*                                                                      *
************************************************************************
                  IF (DoScreen) THEN

                     CALL CWTIME(TCS1,TWS1)

                     ired1 = 1 ! location of the 1st red set
                     Call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,
     &                                  DIAH,DIAG)

                     CALL CWTIME(TCS2,TWS2)
                     tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                     tscrn(2) = tscrn(2) + (TWS2 - TWS1)

                  ENDIF

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      c) 1st MO transformation of DH(mu,nu)                           *
*            Y(mu)[k] = sum_nu  DH(mu,nu) * |C(nu)[k]|                 *
*                                                                      *
************************************************************************

                  nInd = 1
                  Do jDen=1,nKvec
*
** Choose which MO sets on each side
*
                    iMOleft=jDen
                    iMOright=jDen

                    n1 = nIt(iMOright)
                    n2 = nItMx

                    pYik(1:n1,1:n2) => Yik(1:n1*n2)

                    If (DoCAS.and.lSA) iMOright=jDen+2
*

                    Do kSym=1,nSym

                       lSym=MulD2h(JSYM,kSym)
                       Nik= nChOrb_(kSym,iMOleft)*nChOrb_(lSym,iMOright)
                       nIJR(kSym,lSym,jDen) = Nik
                       nIJ1(kSym,lSym,jDen) = Nik
                       If ((JSYM.eq.1).and.iMOleft.eq.iMOright)
     &                                Nik = nChOrb_(kSym,iMOleft)
     &                                    *(nChOrb_(kSym,iMOleft)+1)/2
                       nIJ1(kSym,lSym,jDen) = Nik

                       If (Nik.eq.0) Cycle

                       iS = 1
                       iE = nChOrb_(lSym,iMOright) * JNUM

                       Lik(1:JNUM,1:nChOrb_(lSym,iMOright))=>Aux(iS:iE)

                       iS = iE +1
                       iE = iE + Nik * JNUM

                       Rik(1:Nik*JNUM) => Aux(iS:iE)

                       Rik(:)=Zero

                       Do jK=1,nChOrb_(kSym,iMOleft)
                          jK_a = jK + kOff(kSym,iMOleft)

                           Lik(:,:)=Zero
                           Lab%A0(1:nBas(lSym)*JNUM)=Zero

                        IF (DoScreen .and. iBatch.eq.1) THEN
                           CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------

                           Do ik=1,nBas(kSym)
                              AbsC(ik) = abs(
     &                          MSQ(iMOleft)%SB(kSym)%A2(ik,jK)
     &                                      )
                           End Do

                           If (lSym.ge.kSym) Then

                              mode(1:1)='N'
                              n1 = nBas(lSym)
                              n2 = nBas(kSym)

                           Else ! lSym<kSym

                              mode(1:1)='T'
                              n1 = nBas(kSym)
                              n2 = nBas(lSym)

                           EndIf

                           If (n1>0)
     &                     CALL DGEMV_(Mode(1:1),n1,n2,
     &                                  ONE,DiaH%SB(lSym,kSym)%A2,n1,
     &                                      AbsC,1,
     &                                 ZERO,Ylk(1,jK_a),1)

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      d) 2nd MO transformation of DH(mu,nu)                           *
*            Y(i)[k] = sum_mu  |C(mu)[i]| * Y(mu)[k]                   *
*                                                                      *
************************************************************************

                           If ((kSym.ne.lSym).or.
     &                        (iMOleft.ne.iMOright)) Then
                               iStart=1
                           Else
                               iStart=jK
                           EndIf

                           nQo=0
                           Do i=iStart,nChOrb_(lSym,iMOright)

                              Do ik=1,nBas(lSym)
                                 AbsC(ik) = abs(
     &                             MSQ(iMOright)%SB(lSym)%A2(ik,i)
     &                                         )
                              End Do
*
                              pYik(i,jK_a)=ddot_(nBas(lSym),
     &                                        AbsC,1,Ylk,1)

                              If (pYik(i,jK_a).ge.xtau) Then
                                 nQo=nQo+1
                                 If((iBatch.ne.1) .or.
     &                               (JRED.ne.1)) Go To 1111
                                 nQoT = nQoT + 1
                                 If((lSym .eq. kSym) .and.
     &                              (i .ne. jK)      .and.
     &                              (iMOright.eq.iMOleft)) Then
                                    nQoT = nQoT + 1
                                 End If
                                 If((lSym.eq.kSym).and.
     &                              (iMOleft.eq.iMOright)) Then
                                 End If
 1111                            Indik(1+nQo,nInd)=i
                              Endif

                           End Do
                           Indik(1,nInd)=nQo
************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      e) List the shells present in Y(l)[k] by the largest element    *
*         and sort the list                                            *
*                                                                      *
************************************************************************

                           Do ish=1,nShell
                              YshMax=zero
                              Do ibs=1,nBasSh(lSym,ish)
                                 ibs_a = koffSh(ish,lSym)+ibs
                                 YshMax = Max(YshMax,Ylk(ibs_a,1))
                              End Do
                              MLk(ish) = YshMax
                           End Do

                           Do ish=1,nShell
                              Indx(ish,nInd) = ish
                           End Do

************************************************************************
*                                                                      *
*   1) Screening                                                       *
*                                                                      *
*      f) Screening                                                    *
*                                                                      *
*   Here we use a non-exact bound for the exchange matrix to achieve   *
*   linear scaling. The positive definiteness of the exchange matrix   *
*   combined with the structure of the density matrix makes this       *
*   bound acceptable and likely to be almost exact for what concerns   *
*   the exchange energy                                                *
*                                                                      *
*   The exact bounds (quadratic scaling of the MO transformation)      *
*   would be                                                           *
*      If (MLk(jml)*MLk(1).ge. tau) then                               *
*                                                                      *
*                                                                      *
************************************************************************

                           numSh=0  ! # of significant shells
                           jml=1
                           Do while (jml.le.nShell)

                              YMax=MLk(jml)
                              jmlmax=jml

                              Do iml=jml+1,nShell  ! get the max
                                 If (MLk(iml).gt.YMax) then
                                    YMax = MLk(iml)
                                    jmlmax = iml
                                 Endif
                              End Do

                              If(jmlmax.ne.jml) then  ! swap positions
                                xTmp = MLk(jml)
                                iTmp = Indx(jml,nInd)
                                MLk(jml) = YMax
                                Indx(jml,nInd)=Indx(jmlmax,nInd)
                                MLk(jmlmax) = xTmp
                                Indx(jmlmax,nInd) = iTmp
                              Endif

                              If ( MLk(jml) .ge. xtau ) then
                                numSh = numSh + 1
                              else
                                jml=nShell  ! exit the loop
                              endif

                              jml=jml+1

                           End Do

                           Indx(0,nInd) = numSh

                           CALL CWTIME(TCS2,TWS2)
                           tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                           tscrn(2) = tscrn(2) + (TWS2 - TWS1)
C------------------------------------------------------------------
                        ENDIF    ! Screening setup


                        CALL CWTIME(TCT1,TWT1)

************************************************************************
*                                                                      *
*               E X C H A N G E    T E R M                             *
*                                                                      *
*   2) MO transformation                                               *
*      a) 1st half transformation                                      *
*                                                                      *
*      Transform vectors for shells in the list ML[k]                  *
*                                                                      *
*      Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)   *
*         || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||               *
*                                                                      *
************************************************************************

                        Do iSh=1,Indx(0,nInd)

                           iaSh = Indx(iSh,nInd)

                           Lab%Keep(iaSh,1) = .True.

                           ibcount=0

                           Do ibSh=1,nShell

                              iOffShb = kOffSh(ibSh,kSym)

                              iShp = iTri(iaSh,ibSh)

                              If ( iShp_rs(iShp)<=0) Cycle

                             If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                          nBasSh(lSym,iaSh)*
     &                          nBasSh(kSym,ibSh) .gt. 0
     &                .and. Sqrt(Abs(SumClk%Den(iMOleft)%A2(ibSh,jK_a)*
     &                                     SvShp(iShp_rs(iShp),1) ))
     &                          .ge. thrv )Then

                                ibcount = ibcount + 1

                                IF (lSym.ge.kSym) Then

                                    l1 = 1
                                    If (iaSh<ibSh) l1 = 2

*
**  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
** ------------------------------------
*
                                    Mode(1:1)='N'
                                    n1 = nBasSh(lSym,iaSh)*JNUM
                                    n2 = nBasSh(kSym,ibSh)

                                   Call DGEMV_(Mode(1:1),n1,n2,
     &                    One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,
     &                        MSQ(iMOleft)%SB(kSym)%A2(iOffShb+1:,jK),1,
     &                    ONE,Lab%SB(iaSh,lSym,1)%A,1)

                                Else   ! lSym < kSym

                                   l1 = 1
                                   If (ibSh<iaSh) l1 = 2

*
**  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
** ------------------------------------
*
                                    Mode(1:1)='T'
                                    n1 = nBasSh(kSym,ibSh)
                                    n2 = JNUM*nBasSh(lSym,iaSh)

                                    Call DGEMV_(Mode(1:1),n1,n2,
     &                    One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,
     &                        MSQ(iMOleft)%SB(kSym)%A2(iOffShb+1:,jK),1,
     &                    ONE,Lab%SB(iaSh,lSym,1)%A,1)

                                EndIf

                             EndIf


                           End Do
*
** The following re-assignement is used later on to check if the
** iaSh vector LaJ[k] can be neglected because identically zero
*

                           If (ibcount==0) Lab%Keep(iaSh,1) = .False.

                        End Do

************************************************************************
*                                                                      *
*   2) MO transformation                                               *
*      b) 2nd half transformation                                      *
*                                                                      *
************************************************************************

                        nQo=Indik(1,nInd)

                        Do ir=1,nQo

                          it = Indik(1+ir,nInd)

                          Do iSh=1,Indx(0,nInd)

                             iaSh = Indx(iSh,nInd)

                             If (.NOT.Lab%Keep(iaSh,1)) Cycle

                             iS = kOffsh(iaSh,lSym) + 1
*
                             If (lSym.ge.kSym) Then

**  LJi[k] = sum_a  LaJ[k] * Cai
** ------------------------------
*
                                Mode(1:1)='T'
                                n1 = nBasSh(lSym,iaSh)
                                n2 = JNUM

                             Else   ! lSym < kSym

**   LJi[k] = sum_a  LJa[k] * Cai
** --------------------------------
*
                                Mode(1:1)='N'
                                n1 = JNUM
                                n2 = nBasSh(lSym,iaSh)

                             EndIf

                             CALL DGEMV_(Mode(1:1),n1,n2,
     &                          One,Lab%SB(iaSh,lSym,1)%A,n1,
     &                              MSQ(iMOright)%SB(lSym)%A2(iS:,it),1,
     &                          one,Lik(:,it),1)

                          End Do

*
**  Copy LJi[k] in the standard ordered matrix Lik,J
*
                          If ((jSym.eq.1).and.
     &                       (iMOright.eq.iMOleft)) Then
                             itk = it*(it-1)/2 + jK
                          Else
                             itk = nChOrb_(lSym,iMOright)*(jK-1) + it
                          EndIf
                          call dcopy_(JNUM,Lik(:,it),1,
     &                                     Rik(itk:),Nik)

                        End Do

                        nInd = nInd+1

                        CALL CWTIME(TCT2,TWT2)
                        tmotr2(1) = tmotr(1) + (TCT2 - TCT1)
                        tmotr2(2) = tmotr(2) + (TWT2 - TWT1)


                       End Do  ! loop over k MOs

                       CALL CWTIME(TCT1,TWT1)

************************************************************************
*                                                                      *
*   3) Put to disk                                                     *
*                                                                      *
************************************************************************
                       iAdr = Nik*(JVEC-1)
                       call DDAFILE(LuRVec(lSym,jDen),1,Rik,
     &                                                Nik*JNUM,iAdr)

                       CALL CWTIME(TCT2,TWT2)
                       tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                       tmotr(2) = tmotr(2) + (TWT2 - TWT1)

                    End Do   ! loop over MOs symmetry

                    pYik=>Null()

                  End Do   ! loop over densities

               Call Deallocate_Lab(Lab)
               Call mma_deallocate(Aux)
               Call Deallocate_L_Full(L_Full)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *

C ************  END EXCHANGE CONTRIBUTION  ****************

C --- Diagonals updating. It only makes sense if Nscreen > 0

                  If (Update .and. Nscreen .gt. 0) Then

                     CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
C
C --- subtraction is done in the 1st reduced set
#if defined (_MOLCAS_MPP_)
                     If (Is_Real_Par()) then

                      Do krs=1,nRS
                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) - iiBstR(JSYM,1)
                        Do jvc=1,JNUM
                           DiagJ(jrs) = DiagJ(jrs) + Lrs(krs,jvc)**2
                        End Do
                      End Do

                     Else

                      Do krs=1,nRS
                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set
                        Do jvc=1,JNUM
                           Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                        End Do
                      End Do

                     EndIf

#else
                     Do krs=1,nRS
                        mrs = iiBstR(JSYM,iLoc) + krs
                        jrs = IndRed(mrs,iLoc) ! address in 1st red set
                        Do jvc=1,JNUM
                           Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                        End Do
                     End Do
#endif

                     CALL CWTIME(TCS2,TWS2)
                     tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                     tscrn(2) = tscrn(2) + (TWS2 - TWS1)

                  EndIf

               EndIf ! DoExchange

************************************************************************
************************************************************************
**                                                                    **
**                                                                    **
**    Active term                                                     **
**                                                                    **
**                                                                    **
************************************************************************
************************************************************************
               If (DoCAS) Then

                  CALL CWTIME(TCC1,TWC1)
*
** Set up the skipping flags
** The memory used before for the full-dimension AO-vectors
**     is now re-used to store half and full transformed
**     vectors in the active space
*
                  iSwap = 0  ! Lvb,J are returned
                  Call Allocate_SBA(Laq(1),nAorb,nBas,nVec,JSYM,nSym,
     &                              iSwap)

                  iMO2=1
                  Do iMO1=1,nAdens

*                    iSwap_lxy=5 diagonal blocks are triangular
*                    iSwap_lxy=6 diagonal blocks are square
                     iSwap_lxy=5
                     If (iMO1==2) iSwap_lxy=6
                     Call Allocate_SBA(Lxy,nAorb,nAorb,nVec,JSYM,nSym,
     &                                 iSwap_lxy)


************************************************************************
*                                                                      *
*     MO transformation of Cholesky vectors                            *
*                                                                      *
*         1) Lvb,J = sum_a  C(v,a) * Lab,J                             *
*                                                                      *
************************************************************************

                     kMOs = 1  !
                     nMOs = 1  ! Active MOs (1st set)

                     CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                             JSYM,iSwap,IREDC,nMOs,kMOs,
     &                             Aorb(iMO1),Laq(1),DoRead)

                     if (irc.ne.0) then
                        RETURN
                     endif

************************************************************************
*                                                                      *
*     MO transformation of Cholesky vectors                            *
*                                                                      *
*         2) Lvw,J = sum_b  Lvb,J * C(w,b)                             *
*                                                                      *
************************************************************************
                     If ((JSYM.eq.1).and.(iMO1.eq.iMO2)) Then

                        Do iSymb=1,nSym

                           NAv = nAorb(iSymb)

                           If (NAv<1) Cycle

                           Do JVC=1,JNUM
                             !  triangular blocks
                             CALL DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),
     &                             One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAv,
     &                                      Aorb(iMO2)%SB(iSymb)%A2,NAv,
     &                                Zero,Lxy%SB(iSymb)%A2(:,JVC),NAv)

                          End Do

                        End Do

                     Else

                        Do iSymb=1,nSym

                           iSymv = MulD2h(JSYM,iSymb)
                           NAv = nAorb(iSymv)
                           NAw = nAorb(iSymb) ! iSymb=iSymw

                           If(NAv*NAw.ne.0 .and. iSymv.le.iSymb)Then

                            Do JVC=1,JNUM

                             ! square or rectangular blocks
                             CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                             One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                                      Aorb(iMO2)%SB(iSymb)%A2,NAw,
     &                                 Zero,Lxy%SB(iSymv)%A2(:,JVC),NAv)

                            End Do

                           EndIf

                        End Do

                     EndIf
************************************************************************
*                                                                      *
*     Evaluation of the Z_p_k                                          *
*                                                                      *
*         Z(p){#J} = sum_xy  T(xy,p) * L(xy,{#J})                      *
*                                                                      *
*     T(xy,p) : is stored by compound symmetry JSYM                    *
*               the indices {xy} are stored as PACKED (sym x.le.sym y) *
*                                                                      *
************************************************************************
                     Do iTxy=iMO1,nAdens
                       iAvec=iMO1+iTxy-1
                       Do iSymy=1,nSym

                         iSymx=MulD2h(iSymy,JSYM)

                         If (iSymx.le.iSymy.and.nnA(iSymx,iSymy).ne.0)
     &                      Then

                            ipZp = iOffZp + nnP(JSYM)*(JVEC-1) + 1

                            If (iMO1.eq.iMO2) Then

                              ! diagonal symmetry blocks are triangular
                              CALL DGEMM_('T','N',
     &                           nnP(JSYM),JNUM,nnA(iSymx,iSymy),
     &                       ONE,Txy(ipTxy(iSymx,iSymy,iTxy)),nnP(JSYM),
     &                           Lxy%SB(iSymx)%A2,nnA(iSymx,iSymy),
     &                       ONE,Z_p_k(ipZp,iAvec),nnP(JSYM))

                            Else
*MGD may rearrange the loops

                              Do i=1,nnP(JSYM)
                                 ioff=ipTxy(iSymx,iSymy,iTxy)+
     &                                nnA(iSymx,iSymy)*(i-1)

                                Do j=1,JNUM

*MGD don't work with symmetry
                                  temp=Zero
                                  Do k=0,nAOrb(iSymx)-1
                                    Do l=0,k
                                       temp=temp+0.5d0*
     &                                     Txy(ioff+k*(k+1)/2+l)*
     &                       (Lxy%SB(iSymx)%A2(l+1+nAOrb(iSymx)*k,j)+
     &                        Lxy%SB(iSymx)%A2(k+1+nAOrb(iSymx)*l,j))
                                    End Do
                                  End Do

                                  ij = ipZp -1 + i + nnP(JSYM)*(j-1)

                                  Z_p_k(ij,iAvec)= Z_p_k(ij,iAvec)+temp

                                End Do ! j
                              End Do   ! i

                            EndIf

                         Endif

                       End Do
                     End Do

                     Call Deallocate_SBA(Lxy)
                  End Do

                  CALL CWTIME(TCC2,TWC2)
                  tcasg(1) = tcasg(1) + (TCC2 - TCC1)
                  tcasg(2) = tcasg(2) + (TWC2 - TWC1)

                  Call Deallocate_SBA(Laq(1))


               EndIf  ! DoCAS

************************************************************************
************************************************************************
**                                                                    **
**    Epilogue                                                        **
**                                                                    **
************************************************************************
************************************************************************
*                                                                      *
            END DO  ! end batch loop
*                                                                      *
************************************************************************
*                                                                      *

C --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call mma_deallocate(Drs)
              If(iMp2prpt .eq. 2) Call mma_deallocate(Drs2)
            EndIf

999         Continue

C --- Screening control section
            DoScreen = kscreen.eq.Nscreen

            if (.not.DoScreen) then
                kscreen = kscreen + 1
            else
                kscreen = 1
            endif

            If (DoExchange) Then
#if defined (_MOLCAS_MPP_)
               If (Is_Real_Par() .and. Update .and. DoScreen) Then
                  Call GaDsum(DiagJ,nnBSTR(JSYM,1))
                  Call Daxpy_(nnBSTR(JSYM,1),xone,DiagJ,1,
     &                       Diag(1+iiBstR(JSYM,1)),1)
                  Call Fzero(DiagJ,nnBSTR(JSYM,1))
               EndIf
C--- Need to activate the screening to setup the contributing shell
C--- indices the first time the loop is entered .OR. whenever other nodes
C--- have performed screening in the meanwhile
               If (Is_Real_Par().and..not.DoScreen.and.nVrs.eq.0) Then
                  ntv0=ntv0+1
                  DoScreen = (JRED.lt.myJRED1 .or. ntv0.ge.Nscreen)
                  if (DoScreen) ntv0=0
               EndIf
#endif
            EndIf
*                                                                      *
************************************************************************
*                                                                      *
         END DO   ! loop over red sets
*                                                                      *
************************************************************************
*                                                                      *
         If (DoExchange) Then
           Do jDen=1,nKvec
              Do i=1,nSym
                 Call DACLOS(LuRVec(i,jDen))
              End Do
           End Do
         End If

*                                                                      *
************************************************************************
*                                                                      *
      END DO  ! loop over JSYM
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate a field to be used by Compute_A_jk later
*     since allocations cannot be made at that stage
*                                                                      *
************************************************************************
*                                                                      *
      If(DoExchange) THen
        nIJMax = 0
        Do jDen = 1, nKvec
           Do iSym1 = 1, nSym
              Do iSym2 = 1, nSym
                 nIJMax = max(nIJMax,nIJR(iSym1,iSym2,jDen))
              End Do
           End Do
        End Do
        ljkVec = 2*nIJMax
        Call mma_allocate(VJ,ljkVec,Label='VJ')
      End If

      Call mma_deallocate(iShp_rs)
      If (DoExchange) Then
         Call mma_deallocate(SvShp)
         Call mma_deallocate(kOffSh)
         Call mma_deallocate(Indik)
         Call mma_deallocate(Indx)
         Do i = 1, nDen
            SumClk%Den(i)%A2=>Null()
         End Do
         Call mma_deallocate(SumClk%A0)
         Call mma_deallocate(MLk)
         Call mma_deallocate(Yik)
         Call mma_deallocate(Ylk)
         Call mma_deallocate(AbsC)
         Call Deallocate_NDSBA(DiaH)
#if defined (_MOLCAS_MPP_)
         If (Is_Real_Par().and.Update)CALL mma_deallocate(DiagJ)
#endif
         Call mma_deallocate(Diag)
      EndIf


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1
#ifdef _CD_TIMING_
      ChoGet_CPU = TOTCPU
      ChoGet_Wall = TOTWALL
#endif
*                                                                      *
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky Gradients timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'                                CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB CONTRIB.                 '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'SCREENING OVERHEAD               '
     &                           //'         ',tscrn(1),tscrn(2)
         Write(6,'(2x,A26,2f10.2)')'INACT MO-TRANSFORM VECTORS       '
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'INACT MO-TRANSFORM VECTORS 2     '
     &                           //'         ',tmotr2(1),tmotr2(2)
         Write(6,'(2x,A26,2f10.2)')'ACTIVE CONTRIB.                  '
     &                           //'         ',tcasg(1),tcasg(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      irc  = 0


      Return
      END
