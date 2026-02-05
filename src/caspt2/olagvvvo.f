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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine OLagVVVO(iSym,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,T2AO,
     *                    DIA,DI,FIFA,FIMO,A_PT2,MaxVec_PT2)
      USE iSD_data, only: iSD
      use caspt2_global, only: LuGAMMA,LuCMOPT2,LuAPT2,OLag
      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One
      use caspt2_module, only: IFMSCOUP, IFXMS, IFRMS, IFDW, IFSADREF,
     &                         NSYM, NFRO, NISH, NASH, NSSH,
     &                         NBAS, NBAST, NBMX, NSTATE, JSTATE,
     &                         iRlxRoot
#ifdef _MOLCAS_MPP_
      use caspt2_module, only: NFROT
      use caspt2_global, only: nOLag
      USE Para_Info, ONLY: Is_Real_Par
#endif
C
      implicit none
C
      integer(kind=iwp), intent(in) :: iSym, MaxVec_PT2
      real(kind=wp), intent(in) :: DPT2AO(*), DPT2CAO(*), T2AO(*),
     &  DIA(*), DI(*)
      real(kind=wp), intent(inout) :: FPT2AO(*), FPT2CAO(*), A_PT2(*),
     &  FIFA(*), FIMO(*)
C
      character(len=4096) :: RealName
      logical(kind=iwp) :: DoRys, DoCholesky, is_error, Square
      real(kind=wp), allocatable :: T_hbf(:,:,:,:), vLag(:), WRK1(:),
     &                              WRK2(:)
      integer(kind=iwp), allocatable :: iOffAO(:)
      integer(kind=iwp) :: nBasX(8), KEEP(8), IRC, nSymX, id, lRealName,
     &  iost, iSymI, iSymJ, iSymA, iSymB, nocc, i, nSSDM, nDiff, nSkal,
     &  MaxShlAO, iSh, nBasI, jSh, nBasJ, iBas0, iBas, jBas0, jBas,
     &  iOcc, jOcc, loc1, loc2, j, iRec, nOrbA
      integer(kind=iwp), external :: isFreeUnit
C
C     ----- (VV|VO)
C
C     Compute L_{pq} = (pa|jb) * T_{qj}^{ab}, in particular for
C     (p,q) \in  (virt, inact+act). This operation involves
C     (VV|VO) integrals which are stored in neither disk nor memory.
C     The back-transformation is applied to occupied orbital derivatives
C     of two-electron integrals that have two virtual indices (F, G, H).
C     In principle, the algorithm is to avoid (VV|VO) integrals,
C     i.e. U_{pq} for (p,q) = (vir, inact+act), but can also be applied
C     to U_{pq} for (p,q) = (all, inact+act).
C
      call mma_allocate(vLag,nBasT*nBasT,Label='vLag')
      call mma_allocate(WRK1,nBasT*nBasT,Label='WRK1')
      Call GetOrd(IRC,Square,nSymX,nBasX,KEEP)
C
      vLag(:) = Zero
C
      Call DecideOncholesky(DoCholesky)
      If (DoCholesky) Then
        !! No need to save CMOPT2. Just save A_PT2 and B_PT2.
        !! First, save A_PT2 in LuCMOPT2
        If (IFMSCOUP .and. jState /= 1) Then
          call mma_allocate(WRK2,MaxVec_PT2**2,Label='WRK2')

          ! read A_PT2 from LUAPT2
          id = 0
          call ddafile(LUAPT2, 2, WRK2, MaxVec_PT2**2, id)

          A_PT2(1:MaxVec_PT2**2) = A_PT2(1:MaxVec_PT2**2)
     &      + WRK2(1:MaxVec_PT2**2)
          call mma_deallocate(WRK2)
        End If

        ! For SS-CASPT2 I should write A_PT2 on disk only
        ! for the correct iRlxRoot
        if (jState == iRlxRoot .or. IFMSCOUP) then
          ! write A_PT2 in LUAPT2
          id = 0
          call ddafile(LUAPT2,1,A_PT2,MaxVec_PT2**2,id)
        end if

        !rewind(LuGamma)
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        LuGAMMA = isFreeUnit(LuGAMMA)
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                       'DIRECT','UNFORMATTED',
     &                        iost,.TRUE.,
     &                        nBas(iSym)**2*8,'OLD',is_error)
      End If
      !! 2) Compute ERI (mu rho | nu sigma)
      !! 3) Quarter-transformation of ERI
      !!    (mu rho | j sigma) = sum_{j} C_{nu j} (mu rho | nu sigma)
      !! 4) Contract with AO amplitude
      !!    L_{mu i} = sum_{j,rho sigma} T_{ij}^{mu nu}*(mu rho|j sigma)
      !! the third argument is used wrongly, but should be no problem

      !! D_{mu nu} := DPT2AO
      !! FPT2AO = ((mu nu|rho sigma)-(mu rho|nu sigma)/4
      !!          -(mu sigma|nu rho)/4)*D_{rho sigma}
      isymi = 1
      isymj = 1
      isyma = 1
      isymb = 1
C     nocc = nfro(1)+nish(1)+nash(1)
      nocc = nish(1)+nash(1)
      Call VVVO_Drv(nSym,nBas,nIsh,nFro,KEEP,
     *              iSym,iSymI,iSymA,iSymJ,iSymB,
     *              T2AO,vLag,
     *              nOcc,nBasT,nBMX,
     *              CMOPT2(1+nBasT*nFro(iSymA)),
     *              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *              DIA,DI,FIFA,FIMO)
C
      !! Save the half transformed integrals on disk.
      !! It will be used in drvg1.f etc for gradient.
      If (DoCholesky) Then
        !! Do nothing
        Close (LuGAMMA)
      Else
        !! This is only for conventional calculations!!
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        LuCMOPT2 = isFreeUnit(LuCMOPT2)
        Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
     &                       'DIRECT','UNFORMATTED',
     &                        iost,.FALSE.,1,'OLD',is_error)
        !! First, CMOPT2 has to be saved. The MO coefficient matrix in
        !! grvg1.f may be different from CMOPT2.
        Do i = 1, nBasT*nBasT
          Write (LuCMOPT2) CMOPT2(i)
        End Do
        write (LuCMOPT2) (nIsh(1)+nAsh(1))
        write (LuCMOPT2) (nIsh(2)+nAsh(2))
        write (LuCMOPT2) (nIsh(3)+nAsh(3))
        write (LuCMOPT2) (nIsh(4)+nAsh(4))
        write (LuCMOPT2) (nIsh(5)+nAsh(5))
        write (LuCMOPT2) (nIsh(6)+nAsh(6))
        write (LuCMOPT2) (nIsh(7)+nAsh(7))
        write (LuCMOPT2) (nIsh(8)+nAsh(8))
        write (LuCMOPT2) (nFro(1))
        write (LuCMOPT2) (nFro(2))
        write (LuCMOPT2) (nFro(3))
        write (LuCMOPT2) (nFro(4))
        write (LuCMOPT2) (nFro(5))
        write (LuCMOPT2) (nFro(6))
        write (LuCMOPT2) (nFro(7))
        write (LuCMOPT2) (nFro(8))
        !! Number of state-specific density matrix
        !! It is needed, because the separable contribution to the
        !! electron repulsion integral is different
        nSSDM = 0
        If (.not.IFSADREF) Then
          If (nState == 1) Then
            nSSDM = 0
          Else If (IFDW .or. IFRMS) Then
            !! For (X)DW, use nState density matrix
            nSSDM = nState
          Else If (IFXMS) Then
            !! For XMS, use SA density matrix
            nSSDM = 0
          Else If (IFMSCOUP) Then
            !! For MS, use nState density matrix
            nSSDM = nState
          Else
            !! Otherwise, use one SS density matrix
            nSSDM = 1
          End If
        End If
        write (LuCMOPT2) nSSDM

        Close (LuCMOPT2)

C       write(6,*) "mo saved"
C       call sqprt(CMOPT2,nbast)
C
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        LuGAMMA = isFreeUnit(LuGAMMA)
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                       'DIRECT','UNFORMATTED',
     &                        iost,.TRUE.,
     &                        nOcc*nOcc*8,'OLD',is_error)
        if (is_error) then
        write (u6,*) "Something is wrong in opening LuGamma in olagvvvo"
          call abend
        end if
        !  Setup for shell. Why do I have to call IniSew damn here?
        !  The number of shells should be able to be referred globally.
        nDiff=1
        DoRys=.True.
        Call IniSew(DoRys,nDiff)
        Call Nr_Shells(nSkal)
        Call Setup_iSD()
        !! see Include/info.fh
        Allocate (iOffAO(nSkal+1))
        MaxShlAO = 0
        iOffAO(1) = 0
        Do iSh = 1, nSkal
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          If (nBasI > MaxShlAO) MaxShlAO = nBasI
          iOffAO(iSh+1) = iOffAO(iSh)+nBasI
        End Do
        ! nMax = MaxShlAO*MaxShlAO*nOcc*nOcc

        Allocate (T_hbf(nOcc,nOcc,MaxShlAO,MaxShlAO))
        Do iSh = 1, nSkal
          !! iSD(2,iSh): number of AOs of the shell
          !! iSD(3,iSh): number of cont. func. of the shell
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          Do jSh = 1, nSkal
            nBasJ = iSD(2,jSh)*iSD(3,jSh)
            Do iBas0 = 1, nBasI
              iBas = iOffAO(iSh) + iBas0
              Do jBas0 = 1, nBasJ
                jBas = iOffAO(jSh) + jBas0
                Do iOcc = 1, nOcc
                  Do jOcc = 1, nOcc
                    loc1 = jOcc-1 + (jBas-1)*nOcc
     *                   + (iOcc-1)*nOcc*nBasT
     *                   + (iBas-1)*nOcc*nBasT*nOcc
                    loc2 = iOcc-1 + (iBas-1)*nOcc
     *                   + (jOcc-1)*nOcc*nBasT
     *                   + (jBas-1)*nOcc*nBasT*nOcc
                    T_hbf(iOcc,jOcc,iBas0,jBas0)
     *                = T2AO(loc1+1)+T2AO(loc2+1)
                  End Do
                End Do
                iRec = iBas+nBasT*(jBas-1)
                if (ifmscoup .and. jstate /= 1) then
                  read (lugamma,rec=irec) (wrk1(i),i=1,nocc*nocc)
                  call daxpy_(nocc*nocc,One,wrk2,1,
     *                        t_hbf(1,1,ibas0,jbas0),1)
                end if
                Write (LuGamma,Rec=iRec)
     *            ((T_hbf(i,j,iBas0,jBas0),i=1,nOcc),j=1,nOcc)
              End Do
            End Do
          End Do
        End Do
        Deallocate (iOffAO)
        Deallocate (T_hbf)
        Close (LuGAMMA)
        Call Free_iSD()
        call clssew
      End If
C
      !! 5) L_{ai} = sum_{mu} C_{mu a} L_{mu i}
C     CALL DGEMM_('T','N',nSsh(iSym),nOcc,nBasT,
C    *            One,CMOPT2(1+nBasT*nOcc),nBasT,
C    *                vLag,nBasT,
C    *            One,OLAG(nOCC+1),nOrb(iSymA))
C     write(6,*) "olag before vvvo"
C     call sqprt(olag,nbast)
      nOrbA = nFro(iSymA)+nIsh(iSymA)+nAsh(iSymA)+nSsh(iSymA)
      If (DoCholesky) nOcc = nOrbA-nFro(iSymA)
C     write(6,*) "vLag"
C     call sqprt(vLag,nbast)
      CALL DGEMM_('T','N',nOrbA,nOcc,nBasT,
     *            One,CMOPT2,nBasT,vLag,nBasT,
     *            One,OLAG(nOrbA*nFro(iSymA)+1),nOrbA)
C     write(6,*) "olag after vvvo"
C     call sqprt(olag,nbast)
C
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        if (DoCholesky) call GADSUM(OLag,nOLag)
        If (nFroT == 0) Then
          CALL GADSUM (FPT2AO,nBasT**2)
          CALL GADSUM (FPT2CAO,nBasT**2)
        End If
      End If
#endif
      call mma_deallocate(vLag)
      call mma_deallocate(WRK1)
C
      END SUBROUTINE OLagVVVO
C
C-----------------------------------------------------------------------
C
      !! focktwo_drv.f
      Subroutine VVVO_Drv(nSym,nBas,nAux,nFro,Keep,
     *                    iSym,iSymI,iSymJ,iSymK,iSymL,
     &                    T2AO,vLag,nOcc,nBasT,
     &                    nBMX,CMO,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                    DIA,DI,FIFA,FIMO)
C
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      use Definitions, only: wp, iwp, u6
C
      implicit none

      integer(kind=iwp), intent(in) :: nSym, iSym, iSymI, iSymJ, iSymK,
     &  iSymL, nBas(8), nAux(8), nFro(8), Keep(8), nOcc, nBasT, nBMX
      real(kind=wp), intent(in) :: CMO(*), T2AO(*), DPT2AO(*),
     &  DPT2CAO(*), DIA(*), DI(*)
      real(kind=wp), intent(inout) :: vLag(*), FPT2AO(*), FPT2CAO(*),
     &  FIFA(*), FIMO(*)

      logical(kind=iwp) :: DoCholesky
      real(kind=wp), allocatable :: W1(:), W2(:), WRK(:)

      integer(kind=iwp) :: i, j, LBUF
*
* nAux is the number of occupied orbitals
      DoCholesky=.false.
      Call DecideOnCholesky(DoCholesky)
!     IF (DoCholesky.and.ALGO == 2)THEN
      !! I assume ALGO=2 does not exist (it is not documented and under
      !! debugging, according to rasscf/cho_rasscf_rdinp.f)
!       call abend()
!     end if

      call mma_allocate(W2,NBMX*NBMX,Label='W2')
      call mma_allocate(WRK,nBasT*nBasT,Label='WRK')
*
      Call mma_MaxDBLE(LBUF)
      If (DoCholesky) LBUF = NBMX*NBMX+1
*
* Standard building of the Fock matrix from Two-el integrals
*
      call mma_allocate(W1,LBUF,Label='W1')

      If (LBUF < 1+NBMX**2) Then
         WRITE(u6,*)' FockTwo_Drv Error: Too little memory remains for'
     &     //' the call to FOCKTWO.'
         WRITE(u6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(u6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(u6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
         WRITE(u6,*)'    (All in Real*8-size words)'
         Call  ABEND()
      End If
*
      IF (DoCholesky) Then
        Call VVVOX2(nAux,Keep,iSym,iSymI,iSymJ,iSymK,iSymL,
     *              vLag,CMO,WRK,
     *              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *              FIFA,FIMO)
      Else
        Call VVVOX(nSym,nBas,nFro,Keep,
     *             iSymI,iSymJ,iSymK,iSymL,
     *             T2AO,vLag,CMO,nOcc,nBasT,
     *             LBUF,W1,W2,WRK,
     *             DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *             DIA,DI,FIFA,FIMO)
      End If
      !! vLag must be transposed
      !! In VVVOX(2) subroutines, vLag(p,mu) is constructed.
      WRK(1:nBasT**2) = vLag(1:nBasT**2)
      do i = 1, nbast
        do j = 1, nbast
          vlag(i+nbast*(j-1)) = WRK(j+nbast*(i-1))
        end do
      end do
      call mma_deallocate(WRK)
      call mma_deallocate(W1)
      call mma_deallocate(W2)
*
      End SUBROUTINE VVVO_Drv
C
C-----------------------------------------------------------------------
C
      !! focktwo.f
      SUBROUTINE VVVOX(NSYM,NBAS,NFRO,KEEP,
     *                 iSymI,iSymJ,iSymK,iSymL,
     &                 T2AO,vLag,CMO,nOcc,nBasT,LBUF,X1,X2,WRK,
     *                 DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                 DIA,DI,FIFA,FIMO)

      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One, Half
C
      implicit none
C
      integer(kind=iwp), intent(in) :: NSYM, NBAS(8), NFRO(8), KEEP(8),
     &  iSymI, iSymJ, iSymK, iSymL, nOcc, nBasT, LBUF
      real(kind=wp), intent(in) :: T2AO(nOcc,nBasT,nOcc,nBasT),
     &  CMO(nBasT,*), DPT2AO(*), DPT2CAO(*), DIA(*), DI(*)
      real(kind=wp), intent(inout) :: vLag(nBasT,*), X1(*), X2(*),
     &  WRK(*), FPT2AO(*), FPT2CAO(*), FIFA(*), FIMO(*)

      integer(kind=iwp) ::ISTLT(8), ISTSQ(8), iSym, nB, nB2, nB3, nFroT,
     &  nBasI, KEEPI, nBasJ, KEEPJ, iSymIJ, nBasIJ, nBasK, KEEPK, iSMax,
     &  iSymL_, nBasL, KEEPL, nBasKL, IOPT, LPQ, IPQ, NPQ, IP, JQ, IRC,
     &  ISX, ISF, ISD
C
C     write(6,*) "start vvvox"
C     MUL(I,J)=IEOR(I-1,J-1)+1
      ISTSQ(1)=0
      ISTLT(1)=0
      Do iSym = 2, nSym
        nB  = nBas(iSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(iSym) = ISTSQ(iSym-1) + nB2
        ISTLT(iSym) = ISTLT(iSym-1) + nB3
      End Do
C     write(6,*) "a"
      nFroT = 0
      Do iSym = 1, nSym
        nFroT = nFroT + nFro(iSym)
      End Do
C
      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ == 0) Return
C     write(6,*) "b"

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL_ = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL_ > iSMax) Return !! should not
      nBasL  = nBas(iSymL_)
      KEEPL  = KEEP(iSymL_)
      ! nAuxL  = nAux(iSymL_)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL_) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL == 0) Return
C     write(6,*) "c"
C
      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL /= 0) Return
C     write(6,*) "d"
      !! This will not work when the number of the inactive orbital is 0
C     IF (nAuxI+nAuxJ+nAuxK+nAuxL == 0) Return ! frozen orbitals
C     write(6,*) "e"
C     write(6,*) "nbasij = ", nbasij, 6*13
C     write(6,*) "keep=",keepi,keepj,keepk,keepl
C     write(6,*) "CMO"
C     call sqprt(cmo,nbast)
C
      !! (ij|kl)
C     write(6,*) "doing actual calculation"
      If (iSymI == iSymJ .AND. iSymI == iSymK) Then
C       write(6,*) "aa"
C       write(6,*) "nocc,nbast = ", nocc,nbast
        IOPT=1
        LPQ=0
        IPQ=0
        NPQ=0
C       write(6,*) "nbasi = ", nbasi
        DO IP = 1, nBasI
          DO JQ = 1, IP
C           write(6,*) "ip,jq = ", ip,jq
            IPQ=IPQ+1
            LPQ=LPQ+1
            IF (IPQ > NPQ) THEN
              CALL RDORD(IRC,IOPT,iSymI,iSymJ,iSymK,iSymL,X1,LBUF,NPQ)
              IF(IRC > 1) then
                WRITE(u6,*)' Error return code IRC=',IRC
                WRITE(u6,*)' from RDORD call, in FTWOI.'
                CALL Abend()
              end if
              IOPT=2
              IPQ=1
            ENDIF
C           write(6,*) "do"
            ISX=(IPQ-1)*nBasKL+1
            ISF=ISTLT(iSymI)+LPQ
            ISD=ISTLT(iSymI)+1
            !! coulomb
C           TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD),1)
C           FLT(ISF)=FLT(ISF)+TEMP
            !! exchange
            !! f_{ik} = (ij|kl)*D_{jl}
            !! read (ij|kl) for given ij as (kl)
            !! first  dgemv: (kl)*D_{l,i} -> f_{ki}
            !! second dgemv: (kl)*D_{l,i} -> f_{ki}

            !! For (mu nu), read (mu nu|rho sigma) as I(rho sigma)
            !! (mu j|rho sigma) = sum_{nu} C_{nu j} (mu nu|rho sigma)

            !! L_{mu i} = T_{ij}^{mu nu} * (mu rho | j sigma)
            !!          = T_{ij}^{rho sigma} * C_{nu j}
            !!                 * (mu rho | nu sigma)

            !! IP = mu, JQ = rho
            !! X2 = (nu sigma) --> X2' = C_{nu j}*(nu sigma) = (j sigma)
            !! T_{ij}^{rho sigma} * (j sigma)
            !!   -> U_{i rho} for (mu, rho) pairs

            CALL SQUARE (X1(ISX),X2(1),1,nBasK,nBasL)
C     write(6,*) "ip,jq= ",ip,jq
C     write(6,*) "integral"
C     call sqprt(x2,nbask)
C           If (DoCholesky) Then
C           Else
            !! (mu(ip) rho(jq) | nu sigma) -> (mu(ip) rho(jq) | j sigma)
            call dgemm_('T','N',nOcc,nBasT,nBasT,
     *                  One,CMO,nBasT,X2,nBasT,
     *                  Zero,WRK,nOcc)
C           call dgemm_('T','T',nOcc,nBasT,nBasT,
C    *                  One,CMO,nBasT,X2,nBasT,
C    *                  One,WRK,nOcc)
C           write(6,*) "dgemm finished"
            !! wrk(j,sigma) for the given mu(ip), mu(jq)
            !! T2AO(j,sigma,i,rho) = T_{ij}^{rho sigma}
            !! rather than L_{mu i}, L_{i mu} is computed
            !! L_{i mu} = wrk(j,sigma)*(T2AO(j,sigma,i,rho)
            call dgemv_('t',nOcc*nBasT,nOcc,
     *                  One,t2ao(1,1,1,jq),nOcc*nBasT,wrk,1,
     *                  One,vlag(1,ip),1)
            if (ip /= jq) then
              call dgemv_('t',nOcc*nBasT,nOcc,
     *                    One,t2ao(1,1,1,ip),nOcc*nBasT,wrk,1,
     *                    One,vlag(1,jq),1)
C           call dgemm_('T','T',nOcc,nBasT,nBasT,
C    *                  One,CMO,nBasT,X2,nBasT,
C    *                  One,WRK,nOcc)
C           call dgemv_('n',nOcc*nBasT,nOcc,
C    *                  One,t2ao(1,ip,1,1),nOcc*nBasT,wrk,1,
C    *                  One,vlag(1,jq),1)
            end if
C           End If
C
            !! DPT2AO -> FPT2AO transformation
            !! FPT2 = G(DPT2)
            Call FDGTRF(DPT2AO,FPT2AO)
            Call FDGTRF(DPT2CAO,FPT2CAO)
            If (nFroT /= 0) Then
              Call FDGTRF(DIA,FIFA)
              Call FDGTRF(DI,FIMO)
            End If
            !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2AO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2AO(iSF) = FPT2AO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -Half,X2,nBasK,DPT2AO(iSD),1,
!    *                  One,FPT2AO(iSF),1)
!           If (iP /= jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2AO(iSF) = FPT2AO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -Half,X2,nBasK,DPT2AO(iSD),1,
!    *                    One,FPT2AO(iSF),1)
!           End If
!
!           !! DPT2CAO -> FPT2CAO transformation
!           !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2CAO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -Half,X2,nBasK,DPT2CAO(iSD),1,
!    *                  One,FPT2CAO(iSF),1)
!           If (iP /= jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -Half,X2,nBasK,DPT2CAO(iSD),1,
!    *                    One,FPT2CAO(iSF),1)
!           End If
          End Do
        End Do
      End If

!     ELSE IF ( IS == JS .AND. IS /= KS ) THEN
!       ! CASE 2: Integrals are of symmetry type (II/JJ)
!       ! Coulomb terms need to be accumulated only
!     ELSE IF ( IS == KS .AND. JS == LS ) THEN
!       ! CASE 3: Integrals are of symmetry type (IJ/IJ)
!       ! Exchange terms need to be accumulated only
!         IOPT=1
!         LPQ=0
!         IPQ=0
!         NPQ=0
!         DO IP=1,IB
!           DO JQ=1,JB
!             IPQ=IPQ+1
!             LPQ=LPQ+1
!             IF ( IPQ > NPQ ) THEN
!               CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
!               IF(IRC > 1) then
!                 WRITE(u6,*)' Error return code IRC=',IRC
!                 WRITE(u6,*)' from RDORD call, in FTWOI.'
!                 CALL Abend()
!               end if
!               IOPT=2
!               IPQ=1
!             ENDIF
!             ISX=(IPQ-1)*KLB+1
!             IF ( NFI /= 0 ) THEN
!               ISD=ISTSQ(IS)+(IP-1)*IB+1
!               ISF=ISTSQ(JS)+(JQ-1)*JB+1
!               CALL DGEMV_('N',LB,KB,(-Half*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,One,FSQ(ISF),1)
!             ENDIF
!             IF ( NFJ /= 0 ) THEN
!               ISD=ISTSQ(JS)+(JQ-1)*JB+1
!               ISF=ISTSQ(IS)+(IP-1)*IB+1
!               CALL DGEMV_('T',LB,KB,(-Half*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,One,FSQ(ISF),1)
!             ENDIF
!           End Do
!         End Do
!     End If

      RETURN
C
      Contains
C
      Subroutine FDGTRF(DD,FF)

      implicit none

      real(kind=wp), intent(in) :: DD(*)
      real(kind=wp), intent(inout) :: FF(*)

      real(kind=wp) :: Val
      real(kind=wp), external :: ddot_
      integer(kind=iwp) :: iSF, iSD
C
      !! Coulomb
      Val = DDot_(nBasK*nBasL,X2,1,DD(ISTSQ(iSymI)+1),1)
      iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
      FF(iSF) = FF(iSF) + Val
C
      !! Exchange
      iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
      iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
      CALL DGEMV_('N',nBasK,nBasL,
     *           -Half,X2,nBasK,DD(iSD),1,
     *            One,FF(iSF),1)
      If (iP /= jQ) Then
        iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
        FF(iSF) = FF(iSF) + Val
        iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
        iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
        CALL DGEMV_('N',nBasK,nBasL,
     *             -Half,X2,nBasK,DD(iSD),1,
     *              One,FF(iSF),1)
      End If
C
      End Subroutine FDGTRF
C
      End Subroutine VVVOX
C
C-----------------------------------------------------------------------
C
      !! mabe taken from tracho3.f
      SUBROUTINE VVVOX2(nAux,KEEP,iSym,iSymI,iSymJ,iSymK,iSymL,
     *                  vLag,CMO,WRK,
     *                  DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                  FIFA,FIMO)

      use ChoVec_io, only: NVLOC_CHOBATCH
      use Cholesky, only: InfVec, nDimRS
      use caspt2_global, only: LuGAMMA
      use ChoCASPT2, only: numcho_pt2, NCHSPC, MXNVC
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: IFMSCOUP, NSYM, NFROT, NISH, NASH, NSSH,
     &                         NBAS, NBAST, JSTATE, iRlxRoot, NBTCHES
      use Constants, only: Zero, One, Half, Two

      implicit none

#include "warnings.h"

      integer(kind=iwp), intent(in) :: nAux(8), KEEP(8), iSym, iSymI,
     &  iSymJ, iSymK, iSymL
      real(kind=wp), intent(inout) :: vLag(nBasT,*), WRK(nBasT,nBasT),
     &  FPT2AO(*), FPT2CAO(*), FIFA(*), FIMO(*)
      real(kind=wp), intent(in) :: CMO(nBasT,*), DPT2AO(*), DPT2CAO(*)

      real(kind=wp), allocatable :: CHSPC(:), HTSPC(:), HTVec(:)
      integer(kind=iwp) :: ISTLT(8), ISTSQ(8), iSkip(8), ipWRK(8), jSym,
     &  nAuxT, nB, nB2, nB3, nBasI, KEEPI, nBasJ, KEEPJ, iSymIJ, nBasIJ,
     &  nBasK, KEEPK, iSMax, iSymL_, nBasL, KEEPL, nBasKL, nIshI, nAshI,
     &  nSshI, nOrbI, IBATCH_TOT, JRED1, JRED2, JRED, JSTART, NVECS_RED,
     &  ILOC, IRC, NBATCH, JV1, IBATCH, JNUM, JV2, JREDC, NUMV, MUSED,
     &  iVec, i, j
      real(kind=wp) :: tmp

      Do jSym = 1, nSym
        iSkip(jSym) = 1
        ipWRK(jSym) = 1
      End Do

      ISTSQ(1)=0
      ISTLT(1)=0
      nAuxT = 0
      Do jSym = 2, nSym
        nB  = nBas(jSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(jSym) = ISTSQ(jSym-1) + nB2
        ISTLT(jSym) = ISTLT(jSym-1) + nB3
        nAuxT = nAuxT + nAux(jSym)
      End Do

C     write(6,*) "sym=",isymi,isymj,isymk,isyml
      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
C     write(6,*) "nbasij = ", nbasij
      If (nBasIJ == 0) Return

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL_ = 1+iEor(iSymIJ-1,iSymK-1)
C     write(6,*) "isyml,ismax = ", isyml,ismax
      IF (iSymL_ > iSMax) Return !! should not
      nBasL  = nBas(iSymL_)
      KEEPL  = KEEP(iSymL_)
      ! nAuxL  = nAux(iSymL_)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL) nBasKL = (nBasK*(nBasK+1))/2
C     write(6,*) "nbaskl = ", nbaskl
      If (nBasKL == 0) Return
C
C     write(6,*) "keep=",keepi,keepj,keepk,keepl
      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL /= 0) Return
C     write(6,*) "nAux=",nAuxi,nAuxj,nAuxk,nAuxl
      !! This will not work when the number of the inactive orbital is 0
C     IF (nAuxI+nAuxJ+nAuxK+nAuxL == 0) Return ! frozen orbitals
C
      jSym = iSymJ
      ! kSym = iSymK
      ! lSym = iSymL
      nIshI = nIsh(iSym)
      ! nIshJ = nIsh(jSym)
      ! nIshK = nIsh(kSym)
      ! nIshL = nIsh(lSym)
      nAshI = nAsh(iSym)
      ! nAshJ = nAsh(jSym)
      ! nAshK = nAsh(kSym)
      ! nAshL = nAsh(lSym)
      nSshI = nSsh(iSym)
      ! nSshJ = nSsh(jSym)
      ! nSshK = nSsh(kSym)
      ! nSshL = nSsh(lSym)
      nOrbI = nIshI+nAshI+nSshI
C
C     write(6,*) "nchspc = ", nchspc
      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(HTSPC,NCHSPC,Label='HISPC')
      call mma_allocate(HTVec,nBasT*nBasT,Label='HTVEC')
C
      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym) == 0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,jSym)
      JRED2=InfVec(NumCho_PT2(jSym),2,jSym)
C     write(6,*) "jred1,jred2 = ", jred1,jred2

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED == 0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
* For a reduced set, the structure is known, including
* the mapping between reduced index and basis set pairs.
* The reduced set is divided into suitable batches.
* First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        ! JEND=JSTART+NVECS_RED-1

* Determine batch length for this reduced set.
* Make sure to use the same formula as in the creation of disk
* address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

* Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
C         write(6,*) "ibatch,nbatch = ", ibatch,nbatch
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1
C
C         ----- Construct orbital Lagrangian -----
C
          !! CHSPC       :: (mu nu|P)
          !! HTSPC       :: ( q nu|P)
          !! Bra and Ket :: (ia|P) = T_{ij}^{ab}*(jb|P)
          !! In 1), HT_{i mu,P} = C_{mu a}*(ia|P)
          !!        (ia|P) read from disk
          !! In 3), L_{i mu} = HT_{i nu,P} * (mu nu|P)
          !! Then, L_{pi} = C_{mu p} * L_{i mu}^T
          !! Transpose is done in VVVO_Drv2,
          !! and C_{mu p} is in OLagVVVO
C
          !! 1) Half back-transformation of Bra and Ket density
          !! Read the 3c-2e pseudo-density (in MO), and half transform
          CALL VVVOTRA_RI(CMO,CHSPC,HTSPC,
     *                    JNUM,IBATCH_TOT,IBATCH_TOT,nOrbI)
C
          !! 2) read AO Cholesky vectors,
          !!    then, (strange) reduced form -> squared AO (mu nu|iVec)
          JREDC=JRED
* Read a batch of reduced vectors
          CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)
          IF(NUMV /= JNUM) THEN
            write(u6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(u6,*)' read JNUM vectors. Instead it returned NUMV'
            write(u6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(u6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC /= JRED) THEN
            write(u6,*)' Rats! It was assumed that the Cholesky vectors'
            write(u6,*)' in HALFTRNSF all belonged to a given reduced'
            write(u6,*)' set, but they don''t!'
            write(u6,*)' JRED, JREDC:',JRED,JREDC
            write(u6,*)' Back to the drawing board?'
            write(u6,*)' Let the program continue and see what happens.'
          END IF
C
          !! (strange) reduced form -> squared AO (mu nu|iVec)
          !! is it possible to avoid this transformation?
      ! choptr.fh
          Call R2FIP(CHSPC,WRK,ipWRK(iSym),NUMV,
     *               size(nDimRS),infVec,nDimRS,
     *               nBasT,nSym,iSym,iSkip,irc,JREDC)
C
C           ----- Fock-like transformations (if needed) -----
C
          If (nFroT == 0) Then
            Do iVec = 1, NUMV
              Call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),
     *                    DPT2AO,FPT2AO)
              Call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),
     *                    DPT2CAO,FPT2CAO)
            End Do
          End If
C
          !! 3) Contract with Cholesky vectors
          Call DGemm_('N','T',nOrbI,nBasI,nBasI*JNUM,
     *                Two,HTSPC,nOrbI,CHSPC,nBasI,
     *                One,vLag,nBasI)
C
          !! 4) Construct the 3c-2e pseudo-density in AO
          !! D_{p nu} -> D_{mu nu}
          !! i.e., construct B_PT2, used in ALASKA
          Call DGemm_('N','N',nBasI,nBasI*JNUM,nOrbI,
     *                One,CMO(1,1),nBasI,HTSPC,nOrbI,
     *                Zero,CHSPC,nBasI)
C
          !! 5) Save the 3c-2e pseudo-density in the disk
          !! it may be replaced with ddafile
          Do iVec = 1, NUMV
            If (IFMSCOUP .and. jState /= 1) Then
              Read (LuGamma,Rec=iVec+JV1-1) HTVec(1:nBasI**2)
              Call DaXpY_(nBasI**2,One,
     *                    CHSPC(1+nBasI**2*(iVec-1)),1,
     *                    HTVec,1)
              Write (LuGamma,Rec=iVec+JV1-1) HTVec(1:nBasI**2)
            Else
             if (jState == iRlxRoot .or. IFMSCOUP) then
              Write (LuGamma,Rec=iVec+JV1-1)
     *        CHSPC(1+nBasI**2*(iVec-1):nBasI**2*iVec)
             end if
            End If
          End Do
C
          JV1=JV1+JNUM
        End Do
      End Do
C
      call mma_deallocate(CHSPC)
      call mma_deallocate(HTSPC)
      call mma_deallocate(HTVec)
C
      !! Have to (?) symmetrize Fock-transformed matrices
      If (nFroT == 0) Then
        Do i = 1, nBasI
          Do j = 1, i-1
            tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
            FPT2AO(i+nBasI*(j-1)) = Tmp
            FPT2AO(j+nBasI*(i-1)) = Tmp
            tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
            FPT2CAO(i+nBasI*(j-1)) = Tmp
            FPT2CAO(j+nBasI*(i-1)) = Tmp
            If (nFroT /= 0) Then
              tmp = (FIFA(i+nBasI*(j-1))+FIFA(j+nBasI*(i-1)))*Half
              FIFA(i+nBasI*(j-1)) = Tmp
              FIFA(j+nBasI*(i-1)) = Tmp
              tmp = (FIMO(i+nBasI*(j-1))+FIMO(j+nBasI*(i-1)))*Half
              FIMO(i+nBasI*(j-1)) = Tmp
              FIMO(j+nBasI*(i-1)) = Tmp
            End If
          End Do
        End Do
      End If
C
      Return
C
      Contains
C
      Subroutine FDGTRF(ChoVec,DD,FF)
C
      implicit none
C
      real(kind=wp), intent(in) :: ChoVec(*), DD(*)
      real(kind=wp), intent(inout) :: FF(*)

      real(kind=wp) :: Scal
      real(kind=wp), external :: ddot_
C
      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      Call DaXpY_(nBasI**2,Scal,ChoVec,1,FF,1)
C
      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,
     *            One,ChoVec,nBasI,DD,nBasI,
     *            Zero,WRK,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,
     *           -Half,ChoVec,nBasI,WRK,nBasI,
     *            One,FF,nBasI)
C
      End Subroutine FDGTRF
C
      Subroutine VVVOTRA_RI(CMO,CHSPC_,HTSPC_,NVEC,IBSTA,IBEND,nOrbI)
C
      implicit none
C
      real(kind=wp), intent(in) :: CMO(nBasI,nOrbI)
      !! CHSPC is used as a temporary array
      real(kind=wp), intent(inout) :: CHSPC_(*), HTSPC_(nOrbI,nBasT,*)
      integer(kind=iwp), intent(in) :: NVEC, IBSTA, IBEND, nOrbI

      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) :: IPQ, jVec, nBra
C
      !! BraAI
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,CHSPC_,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nIshI
      Do jVec = 1, NVEC
        ! a. AI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nAshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,
     *                  CMO(1,1+nIshI),nBasI,
     *              Zero,HTSPC_(1,1,jVec),nOrbI)
        ! a. AI -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nIshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,
     *                  CMO(1,1),nBasI,
     *              Zero,HTSPC_(1+nIshI,1,jVec),nOrbI)
      End Do
C
      !! BraSI
      Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,CHSPC_,nBra,
     *                      IBSTA,IBEND)
      IPQ = nIshI*nSshI
      Do jVec = 1, NVEC
        ! b. SI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nSshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,
     *                  CMO(1,1+nIshI+nAshI),nBasI,
     *              One,HTSPC_(1,1,jVec),nOrbI)
        ! a. SI -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nIshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,
     *                  CMO(1,1),nBasI,
     *              Zero,HTSPC_(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do
C
      !! BraSA
      Call Cholesky_Vectors(2,Active,Virtual,JSYM,CHSPC_,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nSshI
      Do jVec = 1, NVEC
        ! d. SA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nSshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,
     *                  CMO(1,1+nIshI+nAshI),nBasI,
     *              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
        ! b. SA -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nAshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,
     *                  CMO(1,1+nIshI),nBasI,
     *              One,HTSPC_(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do
C
      !! BraAA
      Call Cholesky_Vectors(2,Active,Active,JSYM,CHSPC_,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nAshI
      Do jVec = 1, NVEC
        ! b. AA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nAshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,
     *                  CMO(1,1+nIshI),nBasI,
     *              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
        ! c. AA -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nAshI,
     *              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,
     *                  CMO(1,1+nIshI),nBasI,
     *              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
      End Do
C
      End Subroutine VVVOTRA_RI
C
      End Subroutine VVVOX2
C
C-----------------------------------------------------------------------
C
      subroutine getritrfinfo(nnbstr_,maxvec_,n2_)
      use Cholesky, only: infvec_N2, MaxVec, nnBstR
      use definitions, only: wp, iwp
      implicit none
      integer(kind=iwp), intent(out) :: nnbstr_(8,3), maxvec_, n2_
      maxvec_ = maxvec
      n2_     = infvec_n2
      nnbstr_ = nnbstr
      end subroutine getritrfinfo
C
C-----------------------------------------------------------------------
C
      Subroutine R2FIP(CHSPC,WRK,ipWRK,NUMV,l_NDIMRS,INFVEC,
     *                 nDimRS,nBasT,nSym0,iSym,iSkip,irc,JREDC)

      Use Cholesky, only: INFVEC_N2, MaxVec, nnBstR
      use definitions, only: wp, iwp
      use Constants, only: Zero

      implicit none

      real(kind=wp), intent(inout) :: CHSPC(*), WRK(*)
      integer(kind=iwp), intent(in) :: ipWRK(*), NUMV, l_NDIMRS,
     &  INFVEC(MAXVEC,INFVEC_N2,*), nDimRS(nSym0,*), nBasT, nSym0, iSym,
     &  iSkip(8)
      integer(kind=iwp), intent(inout) :: irc, JREDC

      integer(kind=iwp) :: kloc, iVec, lscr, JREDL, ipVecL, jloc
C
C     Transform the reduced form to the full form in place
C
      kloc = 0
      Do iVec = 1, NUMV
        If (l_NDIMRS < 1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        kloc = kloc + lscr
      End Do
C
      ipVecL = 1 + kloc !! lscr*(JNUM-1)
      jloc = (NUMV-1)*nBasT**2+1
      Do iVec = NUMV, 1, -1
        If (l_NDIMRS < 1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        ipVecL = ipVecL - lscr
        WRK(1:nBasT**2) = Zero
        Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,1,
     *                  1,1,1,iSym,JREDC,2,ipWRK,WRK,
     *                  iSkip)
        CHSPC(jloc:jloc+nBasT**2-1) = WRK(1:nBasT**2)
        jloc = jloc-nBasT**2
      End Do
C
      Return
C
      End Subroutine R2FIP
