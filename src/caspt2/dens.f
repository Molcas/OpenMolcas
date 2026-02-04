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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE DENS(IVEC,DMAT,UEFF,U0)
      USE CHOVEC_IO, only: nvloc_chobatch
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: real_shift, imag_shift, sigma_p_epsilon
      use caspt2_global, only: do_grad, do_csf, if_invar, iRoot1,
     *                           iRoot2, if_invaria, if_SSDM,
     *                           CLag,CLagFull,OLag,SLag,
     *                           nCLag,
     *                           DPT2_tot,DPT2C_tot,DPT2_AO_tot,
     *                           DPT2C_AO_tot,DPT2Canti_tot,FIMO_all,
     *                           FIFA_all,OMGDER,jStLag,Weight
      use caspt2_global, only: FIMO, FIFA
      use caspt2_global, only: DREF, DMIX, CMOPT2, TORB, NDREF
      use caspt2_global, only: IDCIEX, IDTCEX
      use PrintLevel, only: debug, verbose
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use EQSOLV, only: IVECX, IVECR, IVECC2
      use ChoCASPT2, only: NumCho_PT2, iALGO, MaxVec_PT2
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: IfChol, IFDENS, IFMSCOUP, IFDW, IFSADREF,
     &                         MAXIT, NSYM, NCONF, NFROT, NISH, NRAS1T,
     &                         NRAS2T, NRAS3T, NASH, NAES, NASHT, NORB,
     &                         NBAS, NBAST, NOSQT, NBSQT, iRlxRoot,
     &                         NSTATE, JSTATE, DENORM, ZETA, ORBIN
      use Constants, only: Zero, One, Two, Half

      implicit none

      integer(kind=iwp), intent(in) :: IVEC
      real(kind=wp), intent(inout) :: DMAT(*)
      real(kind=wp), intent(in) :: UEFF(nState,nState),U0(nState,nState)

      real(kind=wp) :: VECROT(nState)

      real(kind=wp),allocatable :: DPT(:),DSUM(:),DPT2(:),DPT2_AO(:),
     *  DPT2C_AO(:),FPT2(:),FPT2C(:),FPT2_AO(:),FPT2C_AO(:),Trf(:),
     *  WRK1(:),WRK2(:),RDMSA(:,:),RDMEIG(:,:),DEPSA(:,:),
     *  DEPSA_diag(:),DIA(:),DI(:),A_PT2(:),T2AO(:),CLagT(:,:),
     *  EigT(:,:),OMGT(:,:),CI1(:)
      real(kind=wp),allocatable,target :: DPT2Canti_(:),DPT2C(:)
      real(kind=wp),pointer :: DPT2Canti(:)

      integer(kind=iwp) :: NDMAT, NDPT, nDPTAO, ISYM, NO, nAO, IDMOFF,
     &  NI, NA, II, IDM, IT, ITABS, ITTOT, IU, IUTOT, IDRF, IUABS, I, J,
     &  nch, iState, JJ, iStLag, ibk, NumChoTot, nOcc, lT2AO, iSQ, iTR,
     &  nOrbI, ISAV, iBasTr, iBasSq, liBasTr, liBasSq, jBasI, IDSOFF,
     &  IP, IQ, IDSUM, nBasI, iBasI
      real(kind=wp) :: wgt, val, Scal, X
      real(kind=wp) :: CPTF0, CPE, TIOTF0, TIOE, CPTF10, TIOTF10, CPUT,
     &  WALLT

      IF (do_grad) THEN
        !! Set indices for densities and partial derivatives
        Call GradPrep(UEFF,VECROT)
C
C Compute total density matrix as symmetry-blocked array of
C triangular matrices in DMAT. Size of a triangular submatrix is
C  (NORB(ISYM)*(NORB(ISYM)+1))/2.
        NDMAT=0
        NDPT=0
        nDPTAO=0
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          nAO = nBas(iSym)
          NDPT=NDPT+NO**2
          NDMAT=NDMAT+(NO*(NO+1))/2
          nDPTAO = nDPTAO + nAO**2
        END DO
        ! shouldn't be necessary, is already done outside
        DMAT(1:NDMAT) = Zero
C First, put in the reference density matrix.
        IDMOFF=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NORB(ISYM)
          DO II=1,NI
            IDM=IDMOFF+(II*(II+1))/2
            DMAT(IDM)=Two
          END DO
          DO IT=1,NA
            ITABS=NAES(ISYM)+IT
            ITTOT=NI+IT
            DO IU=1,IT
              IUABS=NAES(ISYM)+IU
              IUTOT=NI+IU
              IDRF=(ITABS*(ITABS-1))/2+IUABS
              IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
              DMAT(IDM)=DREF(IDRF)
            END DO
          END DO
           IDMOFF=IDMOFF+(NO*(NO+1))/2
        END DO
*       write(6,*)' DENS. Initial DMAT:'
*       WRITE(*,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
        call mma_allocate(DPT,NDPT,Label='DPT')
        call mma_allocate(DSUM,NDPT,Label='DSUM')
        DPT(:) = Zero
        DSUM(:) = Zero
C
C
C
        !! Modify the solution (T; amplitude), if the real- or
        !! imaginary- shift is utilized. We need both the unmodified (T)
        !! and modified (T+\lambda) amplitudes. \lambda can be obtained
        !! by solving the CASPT2 equation, but it can alternatively
        !! obtained by a direct summation only if CASPT2-D.
        !! iVecX remains unchanged (iVecX = T)
        !! iVecR will be 2\lambda

        !! For MS-CASPT2, calling this subroutine is required.
        !! The lambda-equation is solved without iteration only when
        !! MS-CASPT2-D (shift?). Otherwise, solved iteratively.
        !! After this subroutine, iVecR has multi-state weighted (?)
        !! contributions.
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        Call CASPT2_Res(VECROT)
        CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
        IF (IPRGLB >= verbose) THEN
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(u6,'(a,2f10.2)')" Lambda  : CPU/WALL TIME=", cput,wallt
        END IF
C
C
C
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        !! Diagonal part
        CALL TRDNS2D(iVecX,iVecR,DPT,NDPT,VECROT(JSTATE))
        if (.not.if_invaria) then
          do i = 1, norb(1)
            do j = i+1, norb(1)
              dpt(i+norb(1)*(j-1)) = Zero
              dpt(j+norb(1)*(i-1)) = Zero
            end do
          end do
        end if
        DSUM(:) = DSUM(:) + DPT(:)
*       write(6,*)' DPT after TRDNS2D.'
*       WRITE(*,'(1x,8f16.8)')(dpt(i),i=1,ndpt)
        !! Off-diagonal part, if full-CASPT2
        IF (MAXIT /= 0) THEN
          !! off-diagonal are ignored for CASPT2-D
          DPT(:) = Zero
          CALL TRDNS2O(iVecX,iVecR,DPT,NDPT,VECROT(JSTATE))
          DSUM(:) = DSUM(:) + DPT(:)
        END IF
*       write(6,*)' DPT after TRDNS2O.'
*       WRITE(*,'(1x,8f16.8)')(dpt(i),i=1,ndpt)
        CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
        IF (IPRGLB >= verbose) THEN
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(u6,'(a,2f10.2)')" TRDNS2DO: CPU/WALL TIME=", cput,wallt
        END IF
C
        !! D^PT2 in MO
        call mma_allocate(DPT2,NBSQT,Label='DPT2')
        !! D^PT2(C) in MO
        call mma_allocate(DPT2C,NBSQT,Label='DPT2C')
        !! DPTAO1 (D^PT in AO, but not DPTA-01) couples with
        !! the CASSCF density (assume state-averaged) through ERIs.
        !! This density corresponds to the eigenvalue derivative.
        !! This is sometimes referred to as DPT2(AO) else where.
        call mma_allocate(DPT2_AO,NBSQT,Label='DPT2_AO')
        !! DPTAO2 couples with the inactive density.
        !! This density comes from derivative of the generalized
        !! Fock matrix (see for instance Eq. (24) in the 1990 paper).
        !! This is sometimes referred to as DPT2C(AO) else where.
        call mma_allocate(DPT2C_AO,NBSQT,Label='DPT2C_AO')
        !! DPTAO,DPTCAO,FPTAO,FPTCAO are in a block-squared form
        call mma_allocate(FPT2,NBSQT,Label='FPT2')
        call mma_allocate(FPT2C,NBSQT,Label='FPT2C')
        call mma_allocate(FPT2_AO,NBSQT,Label='FPT2_AO')
        call mma_allocate(FPT2C_AO,NBSQT,Label='FPT2C_AO')
        !! Transformation matrix
        call mma_allocate(Trf,NBSQT,Label='TRFMAT')
        nch=0
        If (IfChol) nch=nvloc_chobatch(1)
        call mma_allocate(WRK1,Max(nBasT**2,nch),Label='WRK1')
        call mma_allocate(WRK2,Max(nBasT**2,nch),Label='WRK2')
        !! state-averaged density
        call mma_allocate(RDMSA,nAshT,nAshT,Label='RDMSA')
        !! Derivative of state-averaged density
        call mma_allocate(RDMEIG,nAshT,nAshT,Label='RDMEIG')
C       write(6,*) "olag before"
C       call sqprt(olag,nbast)
C
        DPT2(:) = Zero
        DPT2C(:) = Zero
        DPT2_AO(:) = Zero
        DPT2C_AO(:) = Zero
        FPT2(:) = Zero
        FPT2C(:) = Zero
        FPT2_AO(:) = Zero
        FPT2C_AO(:) = Zero
        If (.not.IfChol) Then
          FIMO_all(:) = Zero
          FIFA_all(:) = Zero
        End If
        RDMSA(:,:) = Zero
        RDMEIG(:,:) = Zero
C
        CLag(:,:) = Zero
        OLag(:) = Zero
C
        If (nFroT /= 0 .or. .not.if_invaria) Then
          call mma_allocate(DIA,NBSQT,Label='DIA')
          call mma_allocate(DI ,NBSQT,Label='DI')
        Else
          call mma_allocate(DIA,1,Label='DIA')
          call mma_allocate(DI ,1,Label='DI')
        End If
C
        If (do_csf) Then
          call mma_allocate(DPT2Canti_,NBSQT,Label='DPT2Canti')
          DPT2Canti_(:) = Zero
          DPT2Canti => DPT2Canti_
        Else
          DPT2Canti => DPT2C
        End If
C
        !! DPT -> DPT2
        !! Note that DPT2 has the index of frozen orbitals.
        !! Note also that unrelaxed (w/o Z-vector) dipole moments with
        !! frozen orbitals must be wrong.
C       call dcopy_(ndpt,[0.0d+00],0,dpt,1)
        If (nFroT == 0 .and. if_invaria) Then
          DPT2(1:nOsqT) = DSUM(1:nOsqT)
        Else
          Call OLagFro0(DSUM,DPT2)
        End If
C
        !! Construct the transformation matrix
        !! It seems that we have to transform quasi-canonical
        !! to CASSCF orbitals. The forward transformation has been
        !! done in ORBCTL.
        !!   C(PT2) = C(CAS)*X    ->    C(CAS) = C(PT2)*X^T
        !!   -> L(CAS) = X*L(PT2)*X^T
        !! inactive and virtual orbitals are not affected.
        Trf(:) = Zero
        Call CnstTrf(TOrb,Trf)
C       call sqprt(trf,nbast)
C
        !! Construct the density matrix used in the Fock operator
        If (IFSADREF) Then
          WRK1(1:nDRef) = Zero
          Do iState = 1, nState
            wgt = Weight(iState)
            WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*DMix(1:nDRef,iState)
          End Do
        Else
          WRK1(1:nDRef) = DMix(1:nDRef,jState)
        End If
        Call SQUARE(WRK1,RDMSA,1,nAshT,nAshT)
C       write(6,*) "state-averaged density matrix"
C       call sqprt(rdmsa,nasht)
C
C       ----- Construct configuration Lagrangian -----
C
        !! For CI coefficient derivatives (CLag)
        !! Calculate the configuration Lagrangian
        !! This is done in the quasi-canonical basis
        call mma_allocate(DEPSA,nAshT,nAshT,Label='DEPSA')
        DEPSA(:,:) = Zero
        !! Derivative of off-diagonal H0 of <Psi1|H0|Psi1>
        IF (MAXIT /= 0) then
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          Call SIGDER(iVecX,iVecR,VECROT(jState))
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB >= verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(u6,'(a,2f10.2)')
     &        " SIGDER  : CPU/WALL TIME=", cput,wallt
          END IF
        end if
        Call CLagX(1,CLag,DEPSA,VECROT)
C       call test3_dens(clag)
#ifdef _MOLCAS_MPP_
        If (Is_Real_Par()) CALL GADSUM (DEPSA,nAshT**2)
#endif
C       write(6,*) "original depsa"
C       call sqprt(depsa,nasht)
C       write(6,*) "original depsa (sym)"
        do i = 1, nasht
          do j = 1, i-1
            val = (DEPSA(i,j)+DEPSA(j,i))*Half
            DEPSA(i,j) = val
            DEPSA(j,i) = val
          end do
        end do
C       call sqprt(depsa,nasht)
C
        If (NRAS1T+NRAS3T /= 0) Then
          !! The density of the independent pairs (off-diagonal blocks)
          !! should be determined by solving Z-vector, so these blocks
          !! should be removed...?
C         write(6,*) "removing DEPSA of off-diagonal blocks"
C         write(6,*) "before"
C         call sqprt(depsa,nasht)
            Do II = 1, nRAS1T
              Do JJ = nRAS1T+1, nAshT
                DEPSA(II,JJ) = Zero
                DEPSA(JJ,II) = Zero
              End Do
            End Do
            Do II = nRAS1T+1, nRAS1T+nRAS2T
              Do JJ = nRAS1T+nRAS2T+1, nAshT
                DEPSA(II,JJ) = Zero
                DEPSA(JJ,II) = Zero
              End Do
            End Do
C         write(6,*) "after"
C         call sqprt(depsa,nasht)
          IF (IPRGLB >= debug)
     *      write(u6,*) "depsa (sym) after removing off-diagonal blocks"
        Else
          IF (IPRGLB >= debug)
     *      write(u6,*) "depsa (sym)"
        End If
        IF (IPRGLB >= verbose) call sqprt(depsa,nasht)
C
        !! Configuration Lagrangian for MS-CASPT2
        !! This is the partial derivative of the transition reduced
        !! density matrices
        If (IFMSCOUP) Then
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          Call DerHEff(CLag,VECROT)
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB >= verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(u6,'(a,2f10.2)')
     &        " DerHEff : CPU/WALL TIME=", cput,wallt
            write(u6,*)
          END IF
        End If
C
        !! I need to add the derivative of the effective Hamiltonian
        !! for MS-CASPT2, but this is done after orbital Lagrangian.
        !! I just have to have IVECC = T + lambda.
C
        !! If CASPT2 energy is not invariant to rotations in active
        !! orbitals, off-diagonal elements of the density obtained
        !! as DEPSA is incorrect, so remove them. The true density
        !! is computed after everything.
        If (.not.if_invar) Then
          !! But, save the diagonal elements
          call mma_allocate(DEPSA_diag,nAshT,Label='DEPSA_diag')
          Call DCopy_(nAshT,DEPSA,nAshT+1,DEPSA_diag,1)
          !! Clear
          DEPSA(:,:) = Zero
        End If
C       write(6,*) "depsad"
C       call sqprt(depsa,nasht)
C
        !! Transform the quasi-variational amplitude (T+\lambda/2?)
        !! in SR (iVecX) to C (iVecC2)
        !! Note that the contribution is multiplied by two
        !! somewhere else (maybe in olagns?)
        If (real_shift /= Zero .or. imag_shift /= Zero
     &      .OR. sigma_p_epsilon /= Zero .OR. IFMSCOUP) Then
          !! Have to weight the T-amplitude for MS-CASPT2
          IF (IFMSCOUP) THEN
            !! add lambda
            CALL PLCVEC(VECROT(jState),Half,IVECX,IVECR)
            CALL PTRTOC(1,IVECR,IVECC2)
            !! T-amplitude
            Do iStLag = 1, nState
              If (iStLag == jState) Cycle
              Scal = VECROT(iStLag)
              If (ABS(Scal) <= 1.0e-12_wp) Cycle
              Call MS_Res(2,jStLag,iStLag,Scal*Half)
            End Do
            If (do_csf) Then
              !! Prepare for something <\Phi_K^{(1)}|Ers|L>
              ibk = IVECC2
              IVECC2 = 7
              Call RHS_ZERO(IVECC2)
              Do iStLag = 1, nState
                If (iStLag == jState) Cycle
                Scal = UEFF(iStLag,iRoot1)*UEFF(jStLag,iRoot2)
     *               - UEFF(jStLag,iRoot1)*UEFF(iStLag,iRoot2)
                Scal = Scal*Half
                If (ABS(Scal) <= 1.0e-12_wp) Cycle
                Call MS_Res(2,jStLag,iStLag,Scal)
              End Do
              IVECC2 = ibk
            End If
          ELSE
            !! Add lambda to the T-amplitude
            CALL PLCVEC(Half,One,IVECR,IVECX)
            CALL PTRTOC(1,IVECX,IVECC2)
          END IF
        End If
C
C         ipTrfL = 1+nAshT*nBasT+nAshT
C         Call DGemm_('n','N',nAshT,nAshT,nAshT,
C    *                One,Trf(ipTrfL),nBasT,DEPSA,nAshT,
C    *                Zero,dpt2c_ao,nAshT)
C         Call DGemm_('N','t',nAshT,nAshT,nAshT,
C    *                One,dpt2c_ao,nAshT,Trf(ipTrfL),nBasT,
C    *                Zero,DEPSA,nAshT)
C
C       !! Just add DEPSA to DPT2
        Call AddDEPSA(DPT2,DEPSA)
        !! Just transform the density in MO to AO
        CALL DPT2_Trf(DPT,DPT2_AO,CMOPT2,DEPSA,DSUM)
C       call mma_deallocate(DEPSA)
        !! Save the AO density
        !! ... write
C
C       ----- Construct orbital Lagrangian -----
C
        If (nFroT /= 0 .or. .not.if_invaria) Then
          !! If frozen orbitals exist, we need to obtain
          !! electron-repulsion integrals with frozen orbitals to
          !! construct the orbital Lagrangian.
          If (.not.IfChol) Call TRAFRO(1)
C
          !! Get density matrix (DIA) and inactive density
          !! matrix (DI) to compute FIFA and FIMO.
          Call OLagFroD(DIA,DI,RDMSA,Trf)
        End If
C
        !! Construct orbital Lagrangian that comes from the derivative
        !! of ERIs. Also, do the Fock transformation of the DPT2 and
        !! DPT2C densities.
        NumChoTot = 0
        If (IfChol) Then
          Do iSym = 1, nSym
            NumChoTot = NumChoTot + NumCho_PT2(iSym)
          End Do
          call mma_allocate(A_PT2,MaxVec_PT2**2,Label='A_PT2')
          A_PT2(:) = Zero
        Else
          call mma_allocate(A_PT2,1,Label='A_PT2')
        End If
        Do iSym = 1, nSym
          nOcc = nIsh(iSym)+nAsh(iSym)
          If (.not.IfChol .or. iALGO /= 1) Then
            lT2AO = nOcc*nOcc*nBasT*nBasT
            call mma_allocate(T2AO,lT2AO,Label='T2AO')
            T2AO(:) = Zero
          Else
            CALL mma_allocate(T2AO,1,Label='T2AO')
          End If
C
          !! Orbital Lagrangian that comes from the derivative of ERIs.
          !! OLagNS computes only the particle orbitals.
C         write(6,*) "ialgo = ", ialgo
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          If (IfChol .and. iALGO == 1) Then
            CALL OLagNS_RI(iSym,DPT2C,DPT2Canti,A_PT2)
          Else
            CALL OLagNS2(iSym,DPT2C,T2AO)
          End If
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB >= verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(u6,'(a,2f10.2)')
     &        " OLagNS  : CPU/WALL TIME=", cput,wallt
          END IF
C         write(6,*) "DPT2C"
C         call sqprt(dpt2c,nbast)
C
          !! MO -> AO transformations for DPT2 and DPT2C
          If ((.not.IfChol .or. iALGO /= 1)
     *       .or.(nFroT == 0 .and. if_invaria)) Then
            Call OLagTrf(1,iSym,CMOPT2,DPT2,DPT2_AO,WRK1)
            Call OLagTrf(1,iSym,CMOPT2,DPT2C,DPT2C_AO,WRK1)
C           write(6,*) "dpt2"
C           call sqprt(dpt2,nbast)
C           write(6,*) "dpt2ao"
C           call sqprt(dpt2_ao,nbast)
          End If
C
          !! Do some transformations relevant to avoiding (VV|VO)
          !! integrals. Orbital Lagrangian for the hole orbitals are
          !! computed. At the same time, F = G(D) transformations are
          !! also performed for D = DPT2 and DPT2C
          !! The way implemented (what?) is just a shit. I cannot find
          !! FIFA and FIMO for frozen orbitals, so I have to construct
          !! them. Here is the transformation of G(D^inact) and G(D).
          !! FIFA_all and FIMO_all computed in this subroutine
          !! is not yet correct. They are just two-electron after this
          !! subroutine.
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          CALL OLagVVVO(iSym,DPT2_AO,DPT2C_AO,
     *                  FPT2_AO,FPT2C_AO,T2AO,
     *                  DIA,DI,FIFA_all,FIMO_all,
     *                  A_PT2,MaxVec_PT2)
        !   write(6,*) "olag after vvvo"
        !   call sqprt(olag,nbast)
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB >= verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(u6,'(a,2f10.2)')
     &        " OLagVVVO: CPU/WALL TIME=", cput,wallt
          END IF
C     write(6,*) "OLag"
C     do i = 1, 144
C       write(6,'(i3,f20.10)') i,olag(i)
C     end do
C      write(6,*) "fpt2ao"
C     call sqprt(fpt2_ao,12)
C     call abend
C
          !! AO -> MO transformations for FPT2AO and FPT2CAO
          If ((.not.IfChol .or. iALGO /= 1)
     *        .or.(nFroT == 0 .and. if_invaria)) Then
            Call OLagTrf(2,iSym,CMOPT2,FPT2,FPT2_AO,WRK1)
            Call OLagTrf(2,iSym,CMOPT2,FPT2C,FPT2C_AO,WRK1)
          End If
C
          call mma_deallocate(T2AO)
        End Do
        call mma_deallocate(A_PT2)
        !! Add DPTC to DSUM for the correct unrelaxed density
        !! Also, symmetrize DSUM
        Call AddDPTC(DPT2C,DSUM)
C
c
C       write(6,*) "fptao after olagns"
C       call sqprt(fpt2_ao,nbast)
C       write(6,*) "fptcao after olagns"
C       call sqprt(fpt2c_ao,nbast)
C       write(6,*) "olag after olagns"
C       call sqprt(olag,nbast)
C
        !! If frozen orbitals exist, frozen-inactive part of the
        !! unrelaxed PT2 density matrix is computed using the orbital
        !! Lagrangian. Additionally, Fock transformation is also
        !! required.
        If (nFroT /= 0 .or. .not.if_invaria) Then
          !! Compute DPT2 density for frozen-inactive
C         write(6,*) "dpt2 before frozen"
C         call sqprt(dpt2,nbast)
          if (.not.ifchol) then
            !! Construct FIFA and FIMO
            Call OLagFro3(FIFA_all,FIMO_all,WRK1,WRK2)
            !! if possible, canonicalize frozen orbitals, and update
            !! FIMO and Trf
          end if
          !! Save DPT in order to subtract later
          WRK1(1:nDPTAO) = DPT2(1:nDPTAO)
          !! Add explicit FIMO and FIFA contributions. Implicit
          !! contributions are all symmetric in frozen + inactive
          !! orbitals, so they do not contribute to frozen density
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                One,FIMO_all,nBasT,DPT2C,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                One,FIMO_all,nBasT,DPT2C,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                One,FIFA_all,nBasT,DPT2,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                One,FIFA_all,nBasT,DPT2,nBasT,
     *                One,OLAG,nBasT)
C
          !! non-invariant in inactive and secondary
          if(.not.if_invaria) then
            !! Construct the density from orbital Lagrangian
            call caspt2_grad_invaria2(DPT2,OLag)
            !! FIFA contributions from the non-invariant density
            DPT2(1:nDPTAO) = DPT2(1:nDPTAO) - WRK1(1:nDPTAO)
            !! Add the non-invariant contribution to unrelaxed density
            Call AddDPTC(DPT2,DSUM)
            CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                  One,FIFA_all,nBasT,DPT2,nBasT,
     *                  One,OLAG,nBasT)
            CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                  One,FIFA_all,nBasT,DPT2,nBasT,
     *                  One,OLAG,nBasT)
            !! Restore the second-order correlated density
            DPT2(1:nDPTAO) = DPT2(1:nDPTAO) + WRK1(1:nDPTAO)
            DPT2(1:nDPTAO) = WRK1(1:nDPTAO)
          end if
          !! Now, compute pseudo-density using orbital Lagrangian
          !! DSUM does not contain frozen orbitals,
          !! so the properties using this density may be inaccurate
          If(nFroT /= 0) Call OLagFro1(DPT2,OLag)
C
          !! Subtract the orbital Lagrangian added above.
          !! It is computed again in EigDer
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *               -One,FIMO_all,nBasT,DPT2C,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *               -One,FIMO_all,nBasT,DPT2C,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *               -One,FIFA_all,nBasT,WRK1,nBasT,
     *                One,OLAG,nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *               -One,FIFA_all,nBasT,WRK1,nBasT,
     *                One,OLAG,nBasT)
C         write(6,*) "dpt after frozen"
C         call sqprt(dpt2,nbast)
C
          !! Fock transformation for frozen-inactive density
          If (IfChol) Then
            iSym=1
            !! MO -> AO transformations for DPT2 and DPT2C
            Call OLagTrf(1,iSym,CMOPT2,DPT2,DPT2_AO,WRK1)
            Call OLagTrf(1,iSym,CMOPT2,DPT2C,DPT2C_AO,WRK1)
            !! For DF-CASPT2, Fock transformation of DPT2, DPT2C, DIA,
            !! DA is done here, but not OLagVVVO
            !! It seems that it is not possible to do this
            !! transformation in OLagVVVO, because the frozen-part of
            !! the DPT2 is obtained after OLagVVVO.
            FPT2_AO(:) = Zero
            FPT2C_AO(:) = Zero
            Call OLagFro4(1,1,1,1,1,
     *                    DPT2_AO,DPT2C_AO,FPT2_AO,FPT2C_AO,WRK1)
            !! AO -> MO transformations for FPT2AO and FPT2CAO
            Call OLagTrf(2,iSym,CMOPT2,FPT2,FPT2_AO,WRK1)
            Call OLagTrf(2,iSym,CMOPT2,FPT2C,FPT2C_AO,WRK1)
          Else
C           write(6,*) "dpt"
C           call sqprt(dpt2,nbast)
            Call OLagFro2(DPT2,FPT2,WRK1,WRK2)
C         write(6,*) "fpt"
C           call sqprt(dpt2,nbast)
          End If
C         write(6,*) "fifa_all"
C         call sqprt(fifa_all,12)
C         write(6,*) "fimo_all"
C         call sqprt(fimo_all,12)
      !   !! Construct FIFA and FIMO
      !   Call OLagFro3(FIFA_all,FIMO_all,WRK1,WRK2)
C
        Else ! there are no frozen orbitals
          iSQ = 0
          iTR = 0
          Do iSym = 1, nSym
            nOrbI = nOrb(iSym)
            Call SQUARE(FIFA(1+iTR),FIFA_all(1+iSQ),1,nOrbI,nOrbI)
            Call SQUARE(FIMO(1+iTR),FIMO_all(1+iSQ),1,nOrbI,nOrbI)
            iSQ = iSQ + nOrbI*nOrbI
            iTR = iTR + nOrbI*(nOrbI+1)/2
          End Do
        End If
        call mma_deallocate(DIA)
        call mma_deallocate(DI)
C         write(6,*) "fifa_all in dens"
C         call sqprt(fifa_all,nbast)
C         write(6,*) "fimo_all in dens"
C         call sqprt(fimo_all,nbast)
C        write(6,*) "FIFA in natural"
C         Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *                1.0D+00,Trf,nBasT,fifa_all,nBasT,
C    *                0.0D+00,WRK1,nBasT)
C         Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *                1.0D+00,WRK1,nBasT,Trf,nBasT,
C    *                0.0D+00,WRK2,nBasT)
C         call sqprt(wrk2,12)
C
        !! Do some post-process for the contributions that comes from
        !! the above two densities.
C       write(6,*) "olag before eigder"
C       call sqprt(olag,nbast)
C       write(6,*) "fpt2"
C       call sqprt(fpt2,nbast)
        CALL EigDer(DPT2,DPT2C,FPT2_AO,FPT2C_AO,RDMEIG,CMOPT2,
     *              Trf,FPT2,FPT2C,FIFA_all,FIMO_all,RDMSA)
C          call test2_dens(olag,depsa)
C       write(6,*) "olag after eigder"
C       call sqprt(olag,nbast)
C       write(6,*) "Wlag after eigder"
C       call sqprt(wlag,nbast)
C       write(6,*) "rdmeig"
C       call sqprt(rdmeig,nasht)
C       call abend
C
        !! Calculate the configuration Lagrangian again.
        !! The contribution comes from the derivative of eigenvalues.
        !! It seems that TRACI_RPT2 uses CI coefficients of RASSCF,
        !! so canonical -> natural transformation is required.
C       ipTrfL = 1+nAshT*nBasT+nAshT
C       Call DGemm_('N','N',nAshT,nAshT,nAshT,
C    *              1.0D+00,Trf(ipTrfL),nBasT,RDMEIG,nAshT,
C    *              0.0D+00,WRK1,nAshT)
C       Call DGemm_('N','T',nAshT,nAshT,nAshT,
C    *              1.0D+00,WRK1,nAshT,Trf(ipTrfL),nBasT,
C    *              0.0D+00,RDMEIG,nAshT)
        If (.not.if_invar) Then !test
          call mma_allocate(CLagT,nConf,nState,Label='CLagT')
          call mma_allocate(EigT,nAshT,nAshT,Label='EigT')
          CLagT(:,:) = CLag(:,:)
          EigT(:,:) = RDMEIG(:,:)
          If (IFDW .and. zeta >= Zero) then
            call mma_allocate(OMGT,nState,nState,Label='OMGT')
            OMGT(:,:) = OMGDER(:,:)
          end if
        End If
        !! Use canonical CSFs rather than natural CSFs in CLagEig
        ISAV = IDCIEX
        IDCIEX = IDTCEX
        !! Now, compute the configuration Lagrangian
        Call CLagEig(if_SSDM,.false.,CLag,RDMEIG,nAshT)
#ifdef _MOLCAS_MPP_
        If (Is_Real_Par()) CALL GADSUM (CLag,nCLag)
#endif
C
        !! Now, here is the best place to compute the true off-diagonal
        !! active density for non-invariant CASPT2
        If (.not.if_invar) Then
          SLag(:,:) = Zero
          !! Add the density that comes from CI Lagrangian
          Call DEPSAOffC(CLag,DEPSA,FIFA_all,FIMO_all,
     *                   WRK1,WRK2,U0)
          !! Add the density that comes from orbital Lagrangian
          Call DEPSAOffO(OLag,DEPSA,FIFA_all)
          !! Restore the diagonal elements
          Call DCopy_(nAshT,DEPSA_diag,1,DEPSA,nAshT+1)
          call mma_deallocate(DEPSA_diag)
          IF (IPRGLB >= verbose) THEN
            write(u6,*) "DEPSA computed again"
            call sqprt(depsa,nasht)
          END IF
          If (NRAS1T+NRAS3T /= 0) Then
            !! Remove the off-diagonal blocks for RASPT2
            Do II = 1, nRAS1T
              Do JJ = nRAS1T+1, nAshT
                DEPSA(II,JJ) = Zero
                DEPSA(JJ,II) = Zero
              End Do
            End Do
            Do II = nRAS1T+1, nRAS1T+nRAS2T
              Do JJ = nRAS1T+nRAS2T+1, nAshT
                DEPSA(II,JJ) = Zero
                DEPSA(JJ,II) = Zero
              End Do
            End Do
          End If
C         call dcopy_(nasht**2,[0.0d+00],0,depsa,1)
C
          !! We have to do many things again...
          !! Just add DEPSA to DPT2
          Call AddDEPSA(DPT2,DEPSA)
          !! Just transform the density in MO to AO
          CALL DPT2_Trf(DPT,DPT2_AO,CMOPT2,DEPSA,DSUM)
          !! For IPEA shift with state-dependent density
          If (if_SSDM .and. (jState == iRlxRoot .or. IFMSCOUP)) Then
            iSym = 1
            Call OLagTrf(1,iSym,CMOPT2,DPT2,DPT2_AO,WRK1)
          End If
          !! Some transformations similar to EigDer
          Call EigDer2(RDMEIG,Trf,FIFA_all,RDMSA,DEPSA,WRK1,WRK2)
C
          CLag(:,:) = CLagT(:,:) !test
         !test
          RDMEIG(:,:) = RDMEIG(:,:) + EigT(:,:)
          SLag(:,:) = Zero
          call mma_deallocate(CLagT)
          call mma_deallocate(EigT)
          if (IFDW .and. zeta >= Zero) then
            OMGDER(:,:) = OMGT(:,:)
            call mma_deallocate(OMGT)
          end if
          !! RDMEIG contributions
          !! Use canonical CSFs rather than natural CSFs
          !! Now, compute the configuration Lagrangian
          Call CLagEig(if_SSDM,.false.,CLag,RDMEIG,nAshT)
#ifdef _MOLCAS_MPP_
          If (Is_Real_Par()) CALL GADSUM (CLag,nCLag)
#endif
          !! Now, compute the state Lagrangian and do some projections
C         Call CLagFinal(CLag,SLag)
        End If
C
        !! Restore integrals without frozen orbitals, although not sure
        !! this operation is required.
        If ((nFroT /= 0 .or. .not.if_invaria) .and. .not.IfChol)
     *    Call TRAFRO(2)
C
        IDCIEX = ISAV
        !! Canonical -> natural transformation
        IF(ORBIN == 'TRANSFOR') Then
          Do iState = 1, nState
C           Call CLagX_TrfCI(CLag(1+nConf*(iState-1)))
            Call CLagX_TrfCI(CLag(1,iState))
          End Do
        End If
        ! accumulate configuration Lagrangian only for MS,XMS,XDW,RMS,
        ! but not for SS-CASPT2
        if (jState == iRlxRoot .or. IFMSCOUP) then
          CLagFull(1:nConf,1:nState) = CLagFull(1:nConf,1:nState)
     &      + CLag(1:nConf,1:nState)
        end if
C       Call CLagFinal(CLag,SLag)
C
        !! Transformations of DPT2 in quasi-canonical to natural orbital
        !! basis and store the transformed density so that the MCLR
        !! module can use them.
        ! accumulate only if MS,XMS,XDW or RMS calculation
        ! call RecPrt('DPT2 before', '', DPT2_tot, nBast, nBast)
        if (jState == iRlxRoot .or. IFMSCOUP) then
          Call DPT2_TrfStore(One,DPT2,DPT2_tot,Trf,WRK1)
          Call DPT2_TrfStore(Two,DPT2C,DPT2C_tot,Trf,WRK1)
          If (do_csf) Then
            Call DPT2_TrfStore(One,DPT2Canti,DPT2Canti_tot,Trf,WRK1)
          End If
        end if
        ! call RecPrt('DPT2 after', '', DPT2_tot, nBast, nBast)
C       !! Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              One,Trf,nBasT,DPT,nBasT,
C    *              Zero,WRK1,nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              One,WRK1,nBasT,Trf,nBasT,
C    *              Zero,WRK2,nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C         Call DaXpY_(nSQ,One,WRK2(1+iSQ),1,DPT2_tot(1+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C
C       !! Do the same for DPT2C Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              One,Trf,nBasT,DPT2C,nBasT,
C    *              Zero,WRK1,nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              One,WRK1,nBasT,Trf,nBasT,
C    *              Zero,WRK2,nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C        Call DaXpY_(nSQ,Two,WRK2(1+iSQ),1,DPT2C_tot(1+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C       call abend()
C       call sqprt(RDMEIG,nAshT)
C
        !! square -> triangle so that the MCLR module can use the AO
        !! densities. Do this for DPT2AO and DPT2CAO (defined in
        !! caspt2_grad.f and caspt2_grad.h).
        ! accumulate only if MS,XMS,XDW or RMS calculation
        if (jState == iRlxRoot .or. IFMSCOUP) then
          iBasTr = 1
          iBasSq = 1
          Do iSym = 1, nSym
            nBasI = nBas(iSym)
            liBasTr = iBasTr
            liBasSq = iBasSq
            Do iBasI = 1, nBasI
              Do jBasI = 1, iBasI
                liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
                If (iBasI == jBasI) Then
                  DPT2_AO_tot(liBasTr)  = DPT2_AO(liBasSq)
                  DPT2C_AO_tot(liBasTr) = DPT2C_AO(liBasSq)
                Else
                  DPT2_AO_tot(liBasTr) = DPT2_AO(liBasSq)*Two
                  DPT2C_AO_tot(liBasTr) = DPT2C_AO(liBasSq)*Two
                End If
                liBasTr = liBasTr + 1
              End Do
            End Do
            iBasTr = iBasTr + nBasI*(nBasI+1)/2
            iBasSq = iBasSq + nBasI*nBasI
          End Do
        end if
C
        !! If the density matrix used in the Fock operator is different
        !! from the averaged density in the SCF calculation, we need an
        !! additional term for electron-repulsion integral.
        !! Here prepares such densities.
        !! The first one is just DPT2AO, while the second one is the
        !! difference between the SS and SA density matrix. because the
        !! SA density-contribution will be added and should be
        !! subtracted
        ! This should be done only for iRlxRoot
        If (if_SSDM .and. (jState == iRlxRoot .or. IFMSCOUP)) Then
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
!         If (.not.if_invar) Then
!           write(6,*) "SS density matrix with IPEA not implemented"
!           Call abend()
!         End If
C
          !! Construct the SCF density
          !! We first need to construct the density averaged over all
          !! roots involved in SCF.
          WRK1(1:nDRef) = Zero
          call mma_allocate(CI1,nConf,Label='CI1')
          Do iState = 1, nState
            Call LoadCI_XMS('N',1,CI1,iState,U0)
            call POLY1(CI1,nConf)
            call GETDREF(WRK2,nDRef)
            wgt = Weight(iState)
            WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*WRK2(1:nDRef)
          End Do
          call mma_deallocate(CI1)
          !! WRK2 is the SCF density (for nstate=nroots)
          Call SQUARE(WRK1,WRK2,1,nAshT,nAshT)
          Call DaXpY_(nAshT**2,-1.0D+00,WRK2,1,RDMSA,1)
          !! Construct the SS minus SA density matrix in WRK1
          Call OLagFroD(WRK1,WRK2,RDMSA,Trf)
          !! Subtract the inactive part
          WRK1(1:nBasT**2) = WRK1(1:nBasT**2) - WRK2(1:nBasT**2)
          !! Here we should use DPT2_AO??
          !! Save
          If (IfChol) Then
            Call CnstAB_SSDM(DPT2_AO,WRK1)
          Else
            !! Well, it is not working any more. I need to use
            !! Position='APPEND', but it is not possible if I need to
            !! use molcas_open or molcas_open_ext2
            write(u6,*) "It is not possible to perform this calculation"
            write(u6,*) "(non-state averaged density without"
            write(u6,*) "density-fitting or Cholesky decomposition)"
            write(u6,*) "Please use DF or CD"
            call abend()
          End If
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB >= VERBOSE) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(u6,'(a,2f10.2)')
     &        " SSDM    : CPU/WALL TIME=", cput,wallt
          END IF
        End If
C       write(6,*) "pt2ao"
C       call sqprt(DPT2_AO,12)
        call mma_deallocate(DEPSA)
C
        call mma_deallocate(DPT2)
        call mma_deallocate(DPT2C)
        call mma_deallocate(DPT2_AO)
        call mma_deallocate(DPT2C_AO)
        call mma_deallocate(FPT2)
        call mma_deallocate(FPT2C)
        call mma_deallocate(FPT2_AO)
        call mma_deallocate(FPT2C_AO)
        If (do_csf) call mma_deallocate(DPT2Canti_)
        DPT2Canti => null()

        !! Finalize OLag (anti-symmetrize) and construct WLag
        Call OLagFinal(OLag,Trf)
C
        call mma_deallocate(TRF)
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)
        call mma_deallocate(RDMSA)
        call mma_deallocate(RDMEIG)
        DENORM = One
        !! end of with gradient
      ELSE
        !! without gradient
C Compute total density matrix as symmetry-blocked array of
C triangular matrices in DMAT. Size of a triangular submatrix is
C  (NORB(ISYM)*(NORB(ISYM)+1))/2.
        NDMAT=0
        NDPT=0
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          NDPT=NDPT+NO**2
          NDMAT=NDMAT+(NO*(NO+1))/2
        END DO
        DMAT(1:NDMAT) = Zero
C First, put in the reference density matrix.
        IDMOFF=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NORB(ISYM)
          DO II=1,NI
            IDM=IDMOFF+(II*(II+1))/2
            DMAT(IDM)=Two
          END DO
          DO IT=1,NA
            ITABS=NAES(ISYM)+IT
            ITTOT=NI+IT
            DO IU=1,IT
              IUABS=NAES(ISYM)+IU
              IUTOT=NI+IU
              IDRF=(ITABS*(ITABS-1))/2+IUABS
              IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
              DMAT(IDM)=DREF(IDRF)
            END DO
          END DO
           IDMOFF=IDMOFF+(NO*(NO+1))/2
        END DO
*       WRITE(6,*)' DENS. Initial DMAT:'
*       WRITE(6,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
        call mma_allocate(DPT,NDPT,Label='DPT')
        call mma_allocate(DSUM,NDPT,Label='DSUM')
        DSUM(1:NDPT) = Zero

C The 1st order contribution to the density matrix
        DPT(1:NDPT) = Zero
        CALL TRDNS1(IVEC,DPT,NDPT)
        DSUM(1:NDPT) = DSUM(1:NDPT) + DPT(1:NDPT)
*       WRITE(6,*)' DPT after TRDNS1.'
*       WRITE(6,'(1x,8f16.8)')(dpt(i),i=1,ndpt)
        DPT(1:NDPT) = Zero
        CALL TRDNS2D(IVEC,IVEC,DPT,NDPT,1.0D+00)
        IF(IFDENS) THEN
C The exact density matrix evaluation:
          CALL TRDTMP(DPT,NDPT)
        ELSE
C The approximate density matrix evaluation:
          CALL TRDNS2A(IVEC,IVEC,DPT,NDPT)
        END IF
        DSUM(1:NDPT) = DSUM(1:NDPT) + DPT(1:NDPT)
*       WRITE(6,*)' DPT after TRDNS2D.'
*       WRITE(6,'(1x,8f16.8)')(dpt(i),i=1,ndpt)
        DPT(1:NDPT) = Zero
        CALL TRDNS2O(IVEC,IVEC,DPT,NDPT,1.0D+00)
        DSUM(1:NDPT) = DSUM(1:NDPT) + DPT(1:NDPT)
        ! WRITE(6,*)' DPT after TRDNS2O.'
        ! WRITE(6,'(1x,8f16.8)')(dpt(i),i=1,ndpt)
      END IF
C
      call mma_deallocate(DPT)
      IDMOFF=0
      IDSOFF=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        DO IP=1,NO
          DO IQ=1,IP
            IDM=IDMOFF+(IP*(IP-1))/2+IQ
            IDSUM=IDSOFF+IP+NO*(IQ-1)
            DMAT(IDM)=DMAT(IDM)+DSUM(IDSUM)
          END DO
        END DO
        IDMOFF=IDMOFF+(NO*(NO+1))/2
        IDSOFF=IDSOFF+NO**2
      END DO
      call mma_deallocate(DSUM)
C Scale with 1/DENORM to normalize
      X=One/DENORM
      If (do_grad) X=One
      DMAT(1:NDMAT) = X*DMAT(1:NDMAT)

CSVC: For true parallel calculations, replicate the DMAT array
C so that the slaves have the same density matrix as the master.
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        IF (.NOT.KING()) THEN
          DMAT(1:NDMAT) = Zero
        END IF
        CALL GADSUM(DMAT,NDMAT)
      END IF
#endif

      END SUBROUTINE DENS
C
C-----------------------------------------------------------------------
C
      Subroutine CnstTrf(Trf0,Trf)
C
      use caspt2_global, only: TraFro
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NRAS1, NRAS2,
     &                         NRAS3, NASH, NSSH, NDEL, NBAS
      use definitions, only: wp, iwp
      use Constants, only: One
C
      implicit none
C
      real(kind=wp), intent(in) :: Trf0(*)
      real(kind=wp), intent(out) :: Trf(*)

      integer(kind=iwp) :: iSQ, iTOrb, ipTrfL, iSym, nBasI, nFroI,
     &  nIshI, nAshI, nSshI, nDelI, NR1, NR2, NR3, nCor, nVir, I, J,
     &  iIsh, jIsh, IJ, iAsh, jAsh, iSsh, jSsh

      iSQ = 0
      iTOrb = 1 ! LTOrb
      ipTrfL = 0
      Do iSym = 1, nSym
        nBasI = nBas(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nSshI = nSsh(iSym)
        nDelI = nDel(iSym)
        NR1   = nRAS1(iSym)
        NR2   = nRAS2(iSym)
        NR3   = nRAS3(iSym)
        nCor  = nFroI + nIshI
        nVir  = nSshI + nDelI
        ipTrfL = ipTrfL + iSQ
        !! frozen + inactive
C       Do iIsh = 1, nFroI + nIshI
C         Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = 1.0D+00
C       End Do
        !! frozen
        If (IfChol) Then
          Do I = 1, nFroI
            Do J = 1, nFroI
              Trf(ipTrfL+I+nBasI*(J-1))
     *          = TraFro(I+nFroI*(J-1))
            End Do
          End Do
        else
          Do iIsh = 1, nFroI
            Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = One
          End Do
        End If
        !! inactive
        Do I = 1, nIshI
          iIsh = nFroI + I
          Do J = 1, nIshI
            jIsh = nFroI + J
            IJ=I-1+nIshI*(J-1)
            Trf(ipTrfL+iIsh+nBasI*(jIsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nIshI*nIshI
        !! RAS1 space
        Do I = 1, NR1
          iAsh = nCor + I
          Do J = 1, NR1
            jAsh = nCor + J
            IJ=I-1+NR1*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR1*NR1
        !! RAS2 space
        Do I = 1, NR2
          iAsh = nCor + NR1 + I
          Do J = 1, NR2
            jAsh = nCor + NR1 + J
            IJ=I-1+NR2*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR2*NR2
        !! RAS3 space
        Do I = 1, NR3
          iAsh = nCor + NR1 + NR2 + I
          Do J = 1, NR3
            jAsh = nCor + NR1 + NR2 + J
            IJ=I-1+NR3*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR3*NR3
C       call sqprt(trf,12)
      ! !! Active
      ! Do iAsh0 = 1, nAshI
      !   iAsh = nCor + iAsh0
      !   Do jAsh0 = 1, nAshI
      !     jAsh = nCor + jAsh0
C     !     Work(ipTrfL+iAsh-1+nBasI*(jAsh-1))
C    *!       = Work(iTOrb+nIshI*nIshI+iAsh0-1+nAshI*(jAsh0-1))
      !     Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *!       = Trf0(iTOrb+iAsh0-1+nAshI*(jAsh0-1))
      !   End Do
      ! End Do
C       call sqprt(trf,12)
        !! virtual + deleted (deleted is not needed, though)
C       Do iSsh = nOcc+1, nOcc+nVir
C         Trf(ipTrfL+iSsh+nBasI*(iSsh-1)) = 1.0D+00
C       End Do
        Do I = 1, nVir
          iSsh = nCor + nAshI + I
          Do J = 1, nVir
            jSsh = nCor + nAshI + J
            IJ=I-1+nVir*(J-1)
            Trf(ipTrfL+iSsh+nBasI*(jSsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nSshI*nSshI
C       call sqprt(trf,12)
        iSQ = iSQ + nBasI*nBasI

C       n123 = nAshI*nAshI !! just for CAS at present
C       iTOrb = iTOrb + n123 + nSshI*nSshI
C     write(6,*) "transformation matrix"
C     call sqprt(trf,nbasi)
      End Do
C
      Return
C
      End Subroutine CnstTrf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AddDEPSA(DPT2,DEPSA)

      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NASHT, NORB,
     &                         NDEL, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

      real(kind=wp), intent(inout) :: DPT2(*)
      real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, iOrb0,
     &  iOrb2, jOrb0, jOrb2, iOrb, jOrb
      real(kind=wp) :: Val

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2 > 0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nAsh(iSym)
            ! iOrb1 = nIsh(iSym)+iOrb0
            iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              ! jOrb1 = nIsh(iSym)+jOrb0
              jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          = DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          + DEPSA(iOrb0,jOrb0)
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nBas(iSym)-nDel(iSym)
            Do jOrb = 1, iOrb-1
              Val =(DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1))
     *             +DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)))*Half
              DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1)) = Val
              DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine AddDEPSA
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AddDPTC(DPTC,DSUM)

      use caspt2_module, only: NSYM, NFRO, NORB, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

      real(kind=wp), intent(in) :: DPTC(*)
      real(kind=wp), intent(inout) :: DSUM(*)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, iOrb0,
     &  iOrb1, jOrb0, jOrb1, iOrb, jOrb
      real(kind=wp) :: Val

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)!-nDel(iSym)
        If (nOrbI2 > 0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nOrb(iSym)
            iOrb1 = nFro(iSym) + iOrb0
            Do jOrb0 = 1, nOrb(iSym)
              jOrb1 = nFro(iSym) + jOrb0
              DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1))
     *          = DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1))
     *          + DPTC(iMO2+iOrb1-1+nOrbI2*(jOrb1-1))
            End Do
          End Do
          !! Symmetrize DSUM
          Do iOrb = 1, nOrb(iSym)
            Do jOrb = 1, iOrb-1
              Val =(DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1))
     *             +DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)))*Half
              DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1)) = Val
              DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine AddDPTC
C
C-----------------------------------------------------------------------
C
      Subroutine TRAFRO(MODE)

      use caspt2_global, only: CMO, CMO_Internal, CMOPT2, NCMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH, NOSH,
     &                         NSSH, NORB
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: MODE

      integer(kind=iwp) :: nFroTmp(8), nOshTmp(8), nOrbTmp(8)
      integer(kind=iwp) :: jSym

      If (Mode == 1) Then
        Do jSym = 1, nSym
          nFroTmp(jSym) = nFro(jSym)
          nOshTmp(jSym) = nOsh(jSym)
          nOrbTmp(jSym) = nOrb(jSym)
          nOsh(jSym) = nFro(jSym)+nIsh(jSym)+nAsh(jSym)
          nOrb(jSym) = nOsh(jSym)+nSsh(jSym)
          nFro(jSym) = 0
        End Do
      End If
C
      Call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
      CMO=>CMO_Internal
      CMO(:)=CMOPT2(:)
      if (IfChol) then
        call TRACHO3(CMO,NCMO)
      else
        call TRACTL(0)
      end if
      Call mma_deallocate(CMO_Internal)
      nullify(CMO)
C
      If (Mode == 1) Then
        Do jSym = 1, nSym
          nFro(jSym) = nFroTmp(jSym)
          nOsh(jSym) = nOshTmp(jSym)
          nOrb(jSym) = nOrbTmp(jSym)
        End Do
      End If
C
      Return
C
      End Subroutine TRAFRO
C
C-----------------------------------------------------------------------
C
      Subroutine DPT2_TrfStore(Scal,DPT2q,DPT2n,Trf,WRK)

      use caspt2_module, only: NSYM, NORB, NDEL, NBAS
      use Constants, only: Zero, One
      use definitions, only: wp, iwp

      implicit none

      real(kind=wp), intent(in) :: Scal, DPT2q(*), Trf(*)
      real(kind=wp), intent(inout) :: DPT2n(*)
      real(kind=wp), intent(out) :: WRK(*)

      integer(kind=iwp) :: iMO, iSym, nOrbI

      iMO = 1
      Do iSym = 1, nSym
        If (nOrb(iSym) > 0) Then
          nOrbI = nBas(iSym)-nDel(iSym)
          !! Quasi-canonical -> natural transformation of DPT2
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,
     *                One,Trf(iMO),nOrbI,DPT2q(iMO),nOrbI,
     *                Zero,WRK,nOrbI)
          Call DGemm_('N','T',nOrbI,nOrbI,nOrbI,
     *                Scal,WRK,nOrbI,Trf(iMO),nOrbI,
     *                One,DPT2n(iMO),nOrbI)
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do

      Return

      End Subroutine DPT2_TrfStore
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DPT2_Trf(DPT2,DPT2AO,CMO,DEPSA,DSUM)
C
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use Constants, only: Zero, One, Half
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NASHT, NORB,
     &                         NDEL, NBAS, NBSQT
C
      implicit none
C
      real(kind=wp), intent(inout) :: DPT2(*), DPT2AO(*), DSUM(*)
      real(kind=wp), intent(in) :: CMO(*), DEPSA(nAshT,nAshT)

      real(kind=wp), allocatable :: WRK(:)
      real(kind=wp) :: Val
      integer(kind=iwp) :: iCMO, iAO, iMO, iSym, nBasI, nOrbI, iOrb0,
     &  iOrb, jOrb0, jOrb
C
      !! DPT2 transformation
      !! Just transform DPT2 (in MO, block-squared) to DPT2AO (in AO,
      !! block-squared). Also, for DPT2C which couples with the inactive
      !! density matrix.
      call mma_allocate(WRK,NBSQT,Label='WRK')
C
      !! MO -> AO back transformation
      iCMO =1
      iAO = 1
      iMO = 1
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
        If (nORB(ISYM) > 0) Then
          nBasI = nBas(iSym)
          nOrbI = nOrb(iSym)
          !! Add active orbital density
          Do iOrb0 = 1, nAsh(iSym)
            iOrb = nIsh(iSym)+iOrb0
            ! iOrb2= nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              jOrb = nIsh(iSym)+jOrb0
              ! jOrb2= nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + DEPSA(iOrb0,jOrb0)
              DSUM(iMO+iOrb-1+nOrbI*(jOrb-1))
     *          = DSUM(iMO+iOrb-1+nOrbI*(jOrb-1)) + DEPSA(iOrb0,jOrb0)
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nOrb(iSym)
            Do jOrb = 1, iOrb
              Val =(DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *             +DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
            End Do
          End Do
          !! First, DPT2 -> DPT2AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
     *                 One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,
     *                 Zero,WRK,nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
     *                 One,WRK,nBasI,CMO(iCMO),nBasI,
     *                 Zero,DPT2AO(iAO),nBasI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
        iMO  = iMO  + nBasI*nBasI
      End Do
C
      call mma_deallocate(WRK)
C
      END SUBROUTINE DPT2_Trf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer(DPT2,DPT2C,FPT2AO,FPT2CAO,RDMEIG,CMO,Trf,
     *                  FPT2,FPT2C,FIFA,FIMO,RDMSA)
C
      use caspt2_global, only: OLag
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NFROT, NISH, NASH, NASHT,
     &                         NORB, NDEL, NBAS, NBAST, NBSQT
      use Constants, only: Zero, One, Two
C
      implicit none
C
      real(kind=wp), intent(in) :: DPT2(*), DPT2C(*), FPT2AO(*),
     &  FPT2CAO(*), CMO(*), Trf(*), FPT2(*), FPT2C(*), FIFA(*), FIMO(*),
     &  RDMSA(*)
      real(kind=wp), intent(inout) :: RDMEIG(*)

      real(kind=wp),allocatable :: WRK1(:),FPT2_loc(:),FPT2C_loc(:),
     *                             RDMqc(:)

      integer(kind=iwp) :: iCMO, iAO, iSym, nBasI, nOrbI, iSQ, nFroI,
     &  nIshI, nAshI, nCor, iSQA, iT, iTabs, iU, iUabs, iTU, iTUA
C
      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(FPT2_loc,NBSQT,Label='FPT2_loc')
      call mma_allocate(FPT2C_loc,NBSQT,Label='FPT2C_loc')
C
      !! AO -> MO transformation
      iCMO =1
      iAO = 1
      if (nfrot /= 0) then
        FPT2_loc(1:NBSQT) = FPT2(1:NBSQT)
        FPT2C_loc(1:NBSQT) = FPT2C(1:NBSQT)
      else
        DO iSym = 1, nSym
          iCMO = iCMO  + nBas(iSym)*nFro(iSym)
C         iOFF = iWTMP + nBas(iSym)*nBas(iSym)
          If (nOrb(iSym) > 0) Then
            nBasI = nBas(iSym)
            nOrbI = nOrb(iSym)
            !! First, FPT2(AO) -> FPT2(MO)
            CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                   One,CMO(iCMO),nBasI,FPT2AO(iAO),nBasI,
     *                   Zero,WRK1,nOrbI)
            CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                   One,WRK1,nOrbI,CMO(iCMO),nBasI,
     *                   Zero,FPT2_loc(iAO),nOrbI)
            !! Second, FPT2C(AO) -> FPT2C(MO)
            CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                   One,CMO(iCMO),nBasI,FPT2CAO(iAO),nBasI,
     *                   Zero,WRK1,nOrbI)
            CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                   One,WRK1,nOrbI,CMO(iCMO),nBasI,
     *                   Zero,FPT2C_loc(iAO),nOrbI)
          END IF
          iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
          iAO  = iAO  + nBasI*nBasI
C         iMO  = iMO  + nBasI*nBasI
        End Do
      end if
C
      FPT2_loc(:)  = Two*FPT2_loc(:)
      FPT2C_loc(:) = Two*FPT2C_loc(:)
      !! construct Fock in MO
C
      iSQ = 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        ! nSshI = nSsh(iSym)
        ! nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        OLAG(iSQ:iSQ+nOrbI*nCor-1) = OLAG(iSQ:iSQ+nOrbI*nCor-1)
     &    + Two*FPT2_loc(iSQ:iSQ+nOrbI*nCor-1)
        !! Active orbital contributions: (p,q) = (all,act)
        call mma_allocate(RDMqc,nAshI**2,Label='RDMqc')
        !  Construct the active density of the orbital energy
        !  Assume the state-averaged density (SS- and XMS-CASPT2)
C       nSeq = 0
C       Call DCopy_(nAshI*nAshI,[0.0D+00],0,WRK1,1)
C       Do iState = 1, nState
C         Wgt  = DWgt(iState,iState)
C         Wgt  = 1.0D+00/nState
C         Call DaXpY_(nDRef,Wgt,DMIX(:,iState),1,WRK1,1)
C       End Do
        !  RDM of CASSCF
        !  RDMSA is defined by a set of natural orbitals.
        !  Here, we have to transform to a set of quasi-canonical
        !  orbitals (RDMqc), so forward transformation is appropriate.
C       Call SQUARE(WRK1,RDMqc,1,nAshI,nAshI)
        RDMqc(1:nAshI**2) = RDMSA(1:nAshI**2)
        !! nbast?
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              One,Trf(iSQ+nBasT*nCor+nCor),nBasT,RDMqc,nAshT,
     *              Zero,WRK1,nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              One,WRK1,nAshT,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              Zero,RDMqc,nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              One,FPT2_loc(iSQ+nOrbI*nCor),nOrbI,RDMqc,nAshI,
     *              One,OLAG(iSQ+nOrbI*nCor),nOrbI)
        call mma_deallocate(RDMqc)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
C       If (nFroI == 0) Then
C         Call SQUARE(FIFA(iSQ),WRK1,1,nOrbI,nOrbI)
C       Else
C         Call OLagFroSq(iSym,FIFA(iSQ),WRK1)
C       End If
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
!    *              Two,WRK1,nOrbI,DPT2(iSQ),nOrbI,
     *              Two,FIFA(iSQ),nOrbI,DPT2(iSQ),nOrbI,
     *              One,OLAG,nOrbI)
C
        !! explicit derivative of the effective Hamiltonian
        !! dfpq/da = d/da(C_{mu p} C_{nu q} f_{mu nu})
        !!         = f_{mu nu}^a + (C_{mu m} U_{mp} C_{nu q}
        !!                       + C_{mu p} C_{nu m} U_{mq}) * f_{mu nu}
        !!         = f_{mu nu}^a + U_{mp} f_{mq} + U_{mq} f_{pm}
        !! U_{pq}  = f_{pm} df_{qm} + f_{mp} df_{mq}
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
!    *              One,WRK1,nOrbI,DPT2C(iSQ),nOrbI,
     *              One,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              One,OLAG,nOrbI)
        CALL DGEMM_('T','N',nOrbI,nOrbI,nOrbI,
!    *              One,WRK1,nOrbI,DPT2C(iSQ),nOrbI,
     *              One,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              One,OLAG,nOrbI)
C       End If
        !! Implicit derivative of inactive orbitals (DPT2C)
        OLAG(iSQ:iSQ+nOrbI*nCor-1) = OLAG(iSQ:iSQ+nOrbI*nCor-1)
     &    + Two*FPT2C_loc(iSQ:iSQ+nOrbI*nCor-1)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
C
C     ----- CASSCF density derivative contribution in active space
C
      iSQ = 1
      iSQA= 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nCor  = nFroI + nIshI
        Do iT = 1, nAshI
          iTabs = nCor + iT
          Do iU = 1, nAshI
            iUabs = nCor + iU
            iTU = iTabs-1 + nOrbI*(iUabs-1)
            iTUA= iT   -1 + nAshI*(iU   -1)
            RDMEIG(iSQA+iTUA)
     *        = RDMEIG(iSQA+iTUA) + FPT2_loc(iSq+iTU)
C           write(6,'(2i3,f20.10)') it,iu,FPT2(iSq+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C
      call mma_deallocate(WRK1)
      call mma_deallocate(FPT2_loc)
      call mma_deallocate(FPT2C_loc)
C
      END SUBROUTINE EigDer
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer2(RDMEIG,Trf,FIFA,RDMSA,DEPSA,WRK1,WRK2)
C
      use caspt2_global, only: OLag
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NASHT,
     &                         NDEL, NBAS, NBAST, NBSQT
      use Constants, only: Zero, One, Two
C
      implicit none
C
      real(kind=wp), intent(out) :: RDMEIG(*), WRK1(*), WRK2(*)
      real(kind=wp), intent(in) :: Trf(*), FIFA(*), RDMSA(*), DEPSA(*)

      real(kind=wp),allocatable :: FPT2_loc(:),RDMqc(:)
      integer(kind=iwp) :: iSQ, iSym, nOrbI, nFroI, nIshI, nAshI, nCor,
     &  iSQA, iT, iTabs, iU, iUabs, iTU, iTUA
C
      call mma_allocate(FPT2_loc,NBSQT,Label='FPT2_loc')
C
      !! Compute G(D), where D=DEPSA
      Call DEPSATrf(DEPSA,FPT2_loc,WRK1,WRK2)
      FPT2_loc(:) = Two*FPT2_loc(:)
C
      iSQ = 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        ! nSshI = nSsh(iSym)
        ! nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        OLAG(iSQ:iSQ+nOrbI*nCor-1) = OLAG(iSQ:iSQ+nOrbI*nCor-1)
     &    + Two*FPT2_loc(iSQ:iSQ+nOrbI*nCor-1)
        !! Active orbital contributions: (p,q) = (all,act)
        call mma_allocate(RDMqc,nAshI**2,Label='RDMqc')
        RDMqc(1:nAshT**2) = RDMSA(1:nAshT**2)
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              One,Trf(iSQ+nBasT*nCor+nCor),nBasT,RDMqc,nAshT,
     *              Zero,WRK1,nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              One,WRK1,nAshT,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              Zero,RDMqc,nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              One,FPT2_loc(iSQ+nOrbI*nCor),nOrbI,RDMqc,nAshI,
     *              One,OLAG(iSQ+nOrbI*nCor),nOrbI)
        call mma_deallocate(RDMqc)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
        CALL DGEMM_('N','T',nOrbI,nAshI,nAshI,
     *              Two,FIFA(1+nOrbI*nCor),nOrbI,DEPSA,nAshI,
     *              One,OLAG(iSQ+nOrbI*nCor),nOrbI)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
C
C     ----- CASSCF density derivative contribution in active space
C
      iSQ = 1
      iSQA= 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nCor  = nFroI + nIshI
        Do iT = 1, nAshI
          iTabs = nCor + iT
          Do iU = 1, nAshI
            iUabs = nCor + iU
            iTU = iTabs-1 + nOrbI*(iUabs-1)
            iTUA= iT   -1 + nAshI*(iU   -1)
            RDMEIG(iSQA+iTUA) = FPT2_loc(iSq+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C
      call mma_deallocate(FPT2_loc)
C
      END SUBROUTINE EigDer2
C
C-----------------------------------------------------------------------
C
      Subroutine DEPSATrf(DEPSA,FPT2,WRK1,WRK2)
C
      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH, NASHT,
     &                         NBAS, NBAST, NBSQT
      use Constants, only: Zero, One, Half
C
      implicit none
C
      real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)
      real(kind=wp), intent(out) :: FPT2(*), WRK1(*), WRK2(*)

      real(kind=wp),allocatable :: DAO(:),DMO(:)
      integer(kind=iwp) :: iSym, iSymA, iSymI, iSymB, iSymJ,
     &  nCorI, nBasI, iAsh, jAsh, iAshI, iOrb, jAshI, jOrb
C
      FPT2(1:nBasT**2) = Zero
C
      iSym = 1
      iSymA= 1
      iSymI= 1
      iSymB= 1
      iSymJ= 1
C
C     If (nFroT.ne.0.and.IfChol) Then
      If (IfChol) Then
        !! DEPSA(MO) -> DEPSA(AO) -> G(D) in AO -> G(D) in MO
        !! The Cholesky vectors do not contain frozen orbitals...
        call mma_allocate(DAO,NBSQT,Label='DAO')
        call mma_allocate(DMO,NBSQT,Label='DMO')
        !! First, MO-> AO transformation of DEPSA
        Do iSym = 1, nSym
          DMO(:) = Zero
          nCorI = nFro(iSym)+nIsh(iSym)
          nBasI = nBas(iSym)
          Do iAsh = 1, nAsh(iSym)
            Do jAsh = 1, nAsh(iSym)
              DMO(nCorI+iAsh+nBasI*(nCorI+jAsh-1)) = DEPSA(iAsh,jAsh)
            End Do
          End Do
          Call OLagTrf(1,iSym,CMOPT2,DMO,DAO,WRK1)
        End Do
        !! Compute G(D)
        WRK1(1:NBSQT) = Zero
        DMO(:) = Zero
        !! it's very inefficient
        Call OLagFro4(1,1,1,1,1,
     *                DAO,WRK1,DMO,WRK1,WRK2)
        !! G(D) in AO -> G(D) in MO
        Do iSym = 1, nSym
          Call OLagTrf(2,iSym,CMOPT2,FPT2,DMO,WRK1)
        End Do
        call mma_deallocate(DAO)
        call mma_deallocate(DMO)
      Else
        nCorI = nFro(iSym)+nIsh(iSym)
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            FPT2(1:nBasT**2) = FPT2(1:nBasT**2)
     &        + DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)
C
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            FPT2(1:nBasT**2) = FPT2(1:nBasT**2)
     &        - Half*DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)
          End Do
        End Do
      End If
C
      Return
C
      End Subroutine DEPSATrf
C
C-----------------------------------------------------------------------
C
      Subroutine CnstAB_SSDM(DPT2AO,SSDM)
C
      use ChoVec_io, only: NVLOC_CHOBATCH
      use Cholesky, only: InfVec, nDimRS
      use caspt2_global, only: LuGAMMA,LuAPT2
      use ChoCASPT2, only: NumCho_PT2, MaxVec_PT2, NCHSPC, MXNVC
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: NSYM, NBAS, NBAST, NBSQT, NBTCHES
      use Constants, only: Zero, One, Two
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, nProcs, myRank
#endif
C
      implicit none
C
#include "warnings.h"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
C
      real(kind=wp), intent(in) :: DPT2AO(*),SSDM(*)

      integer(kind=iwp) :: iSkip(8), ipWRK(8), nnbstr(8,3)
      character(len=4096) :: RealName
      logical(kind=iwp) :: is_error
      real(kind=wp),allocatable :: A_PT2(:),CHSPC(:),HTVec(:),WRK(:),
     *  V1(:),V2(:),B_SSDM(:)

      integer(kind=iwp) :: maxvec, n2, iSym, NumChoTot, jSym, nBasI, id,
     &  lRealName, iost, JRED1, JRED2, i, ipV1, ipV2, IBATCH_TOT, JRED,
     &  JSTART, NVECS_RED, ILOC, IRC, NBATCH, JV1, IBATCH, JNUM,
     &  JV2, JREDC, NUMV, ipVecL, iVec, lscr, JREDL, MUSED, NUMVI, KV1,
     &  JBATCH_TOT, JBATCH, KNUM, KV2, NUMVJ
      real(kind=wp), external :: ddot_
      integer(kind=iwp), external :: isFreeUnit

#ifdef _MOLCAS_MPP_
      logical(kind=iwp) :: bStat
      integer(kind=iwp), allocatable :: map2(:)
      integer(kind=iwp) :: lg_V1, ILOV1, IHIV1, JLOV1, JHIV1, MV1,
     &  NDIM1, JVG, iRank
#endif
C
      ! INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      call getritrfinfo(nnbstr,maxvec,n2)
      iSym = 1 !! iSym0
      i = 0
C
      NumChoTot = 0
      Do jSym = 1, nSym
        NumChoTot = NumChoTot + NumCho_PT2(jSym)
      End Do
      Do jSym = 1, nSym
        iSkip(jSym) = 1
#ifdef _MOLCAS_MPP_
        !! I do not know why this is necessary
        If (is_real_par()) then
          ipWRK(jSym) = 1
        else
#endif
          ipWRK(jSym) = 1
#ifdef _MOLCAS_MPP_
        end if
#endif
      End Do
C
      nBasI  = nBas(iSym)
C
      call mma_allocate(A_PT2,MaxVec_PT2**2,Label='A_PT2')

      ! Read A_PT2 from LUAPT2
      id = 0
      call ddafile(LUAPT2, 2, A_PT2, MaxVec_PT2**2, id)

      !! Open B_PT2
      Call PrgmTranslate('GAMMA',RealName,lRealName)
      LuGAMMA = isFreeUnit(LuGAMMA)
      Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                     'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      nBas(iSym)**2*8,'OLD',is_error)
C
      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(HTVec,nBasT*nBasT,Label='HTVec')
      call mma_allocate(WRK,nBasT**2,Label='WRK')
      !! V(P) = (mu nu|P)*D_{mu nu}
      call mma_allocate(V1,MaxVec_PT2,Label='V1')
      call mma_allocate(V2,MaxVec_PT2,Label='V2')
      !! B_SSDM(mu,nu,P) = D_{mu rho}*D_{nu sigma}*(rho sigma|P)
      !! Add one more vector
      call mma_allocate(B_SSDM,NCHSPC+NBSQT,Label='B_SSDM')

      IF(NUMCHO_PT2(iSym) == 0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,iSym)
      JRED2=InfVec(NumCho_PT2(iSym),2,iSym)
C
#ifdef _MOLCAS_MPP_
      If (is_real_par()) then
        call mma_allocate(MAP2,nProcs,Label='MAP2')
        MAP2(:) = 0
        MAP2(myRank+1) = NumChoTot ! MJRED2-JRED1+1
        call GAIGOP(MAP2,NPROCS,'+')
C       ndim2 = sum(map2)

        do i = nprocs, 2, -1
          map2(i) = sum(map2(1:i-1))+1
        end do
        map2(1) = 1
        ipV1 = map2(myRank+1)
        ipV2 = map2(myRank+1)
        V1(:) = Zero
        V2(:) = Zero

        bStat = GA_CREATE_IRREG(MT_DBL,nBasT**2,MaxVec_PT2,'WRK',
     *                          1,1,MAP2,NPROCS,lg_V1)
        CALL GA_DISTRIBUTION(LG_V1,MYRANK,ILOV1,IHIV1,JLOV1,JHIV1)
        CALL GA_ACCESS(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,MV1,NDIM1)
      else
#endif
        ipV1 = 1
        ipV2 = 1
#ifdef _MOLCAS_MPP_
      end if
#endif
C
      !! Prepare density matrix
      !! subtract the state-averaged density matrix
C
      IBATCH_TOT=NBTCHES(iSym)
C
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
          ipVecL = 1
          Do iVec = 1, NUMV
C
            !! reduced form -> squared AO vector (mu nu|iVec)
            ! If (l_NDIMRS < 1) Then
            If (size(nDimRS) < 1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              ! lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
              lscr  = nDimRS(iSym,JREDL)
            End If
            WRK(:) = Zero
            Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,1,
     *                      1,1,1,iSym,JREDC,2,ipWRK,WRK,
     *                      iSkip)
            ipVecL = ipVecL + lscr
C
            V1(ipV1+JV1+iVec-2) = DDot_(nBasI**2,DPT2AO,1,WRK,1)
            V2(ipV2+JV1+iVec-2) = DDot_(nBasI**2,SSDM  ,1,WRK,1)
C
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  One,DPT2AO,nBasI,WRK,nBasI,
     *                  Zero,HTVec,nBasI)
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  One,HTVec,nBasI,SSDM,nBasI,
     *                  Zero,B_SSDM(1+nBasT**2*(iVec-1)),nBasI)
          End Do
          NUMVI = NUMV
C
#ifdef _MOLCAS_MPP_
          If (is_real_par()) then
            Call DCopy_(nBasT**2*NUMV,B_SSDM,1,
     *                  DBL_MB(mV1+NDIM1*(JV1-1)),1)
          else
#endif
            KV1=JSTART
            JBATCH_TOT=NBTCHES(iSym)
            DO JBATCH=1,NBATCH
              JBATCH_TOT=JBATCH_TOT+1

              KNUM=NVLOC_CHOBATCH(JBATCH_TOT)
              KV2=KV1+KNUM-1

              JREDC=JRED
              CALL CHO_VECRD(CHSPC,NCHSPC,KV1,KV2,iSym,
     &                              NUMV,JREDC,MUSED)
              Call R2FIP(CHSPC,WRK,ipWRK,NUMV,
     *                   size(nDimRS),infVec,nDimRS,
     *                   nBasT,nSym,iSym,iSkip,irc,JREDC)
C
              !! Exchange part of A_PT2
              NUMVJ = NUMV
              CALL DGEMM_('T','N',NUMVI,NUMVJ,nBasT**2,
     *                   -One,B_SSDM,nBasT**2,CHSPC,nBasT**2,
     *                    One,A_PT2(JV1+MaxVec_PT2*(KV1-1)),
     *                        MaxVec_PT2)
              KV1=KV1+KNUM
            END DO
#ifdef _MOLCAS_MPP_
          end if
#endif
C
          !! Read, add, and save the B_PT2 contribution
          Do iVec = 1, NUMVI
            Read  (Unit=LuGAMMA,Rec=JV1+iVec-1) WRK(1:nBasT**2)
            !! The contributions are doubled,
            !! because halved in PGet1_RI3?
            !! Coulomb
            WRK(1:nBasT**2) = WRK(1:nBasT**2)
     &        + V2(ipV2+iVec-1)*DPT2AO(1:nBasT**2)
     &        + V1(ipV1+iVec-1)*SSDM(1:nBasT**2)
            !! Exchange
            WRK(1:nBasT**2) = WRK(1:nBasT**2)
     &        - B_SSDM(1+nBasT**2*(iVec-1):1+nBasT**2*iVec)
            Write (Unit=LuGAMMA,Rec=JV1+iVec-1) WRK(1:nBasT**2)
          End Do
          JV1=JV1+JNUM
        End Do
      End Do
C
#ifdef _MOLCAS_MPP_
      If (is_real_par()) then
        !! Parallel for the exchange part of A_PT2
        A_PT2(:) = A_PT2/real(NPROCS,kind=wp)
        IBATCH_TOT=NBTCHES(iSym)
        DO JRED=JRED1,JRED2
          CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
          IF(NVECS_RED == 0) Cycle
          ILOC=3
          CALL CHO_X_SETRED(IRC,ILOC,JRED)
          NBATCH=1+(NVECS_RED-1)/MXNVC
* Loop over IBATCH
          JV1=JSTART
          DO IBATCH=1,NBATCH
C           write(6,*) "ibatch,nbatch = ", ibatch,nbatch
            IBATCH_TOT=IBATCH_TOT+1

            JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
            JV2=JV1+JNUM-1

            JREDC=JRED
* Read a batch of reduced vectors
            CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,NUMV,JREDC,MUSED)
            Call R2FIP(CHSPC,WRK,ipWRK,NUMV,size(nDimRS),infVec,nDimRS,
     *                 nBasT,nSym,iSym,iSkip,irc,JREDC)
            NUMVJ = NUMV
C
            !! Exchange part of A_PT2
            JVG   = JV1+MAP2(myRank+1)-1
            Do iRank = 0, NPROCS-1
              CALL GA_DISTRIBUTION(LG_V1,iRank,ILOV1,IHIV1,JLOV1,JHIV1)
              CALL GA_GET(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,B_SSDM,
     *                    NDIM1)
              NUMVI = JHIV1-JLOV1+1
              CALL DGEMM_('T','N',NUMVI,NUMVJ,nBasT**2,
     *                   -One,B_SSDM,nBasT**2,CHSPC,nBasT**2,
     *                    One,A_PT2(JLOV1+MaxVec_PT2*(JVG-1)),
     *                        MaxVec_PT2)
            end do
            JV1=JV1+JNUM
          End Do
        End Do
        CALL GADSUM (A_PT2,MaxVec_PT2**2)
        bStat = GA_Destroy(lg_V1)
C
        CALL GADSUM (V1,MaxVec_PT2)
        CALL GADSUM (V2,MaxVec_PT2)
        call mma_deallocate(MAP2)
      end if
#endif
C
      !! Coulomb for A_PT2
      !! Consider using DGER?
      Call DGEMM_('N','T',MaxVec_PT2,MaxVec_PT2,1,
     *            Two,V1,MaxVec_PT2,V2,MaxVec_PT2,
     *            One,A_PT2,MaxVec_PT2)
C
      ! write to A_PT2 in LUAPT2
      id = 0
      call ddafile(LUAPT2, 1, A_PT2, MaxVec_PT2**2, id)
C
      !! close B_PT2
      Close (LuGAMMA)

      call mma_deallocate(A_PT2)
C
      call mma_deallocate(CHSPC)
      call mma_deallocate(HTVec)
      call mma_deallocate(WRK)
      call mma_deallocate(V1)
      call mma_deallocate(V2)
      call mma_deallocate(B_SSDM)
C
      End Subroutine CnstAB_SSDM
