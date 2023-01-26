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
      SUBROUTINE DENS(IVEC,DMAT,UEFF)
      USE CHOVEC_IO
      use caspt2_output, only: iPrGlb, verbose, debug
      use caspt2_global, only: real_shift, imag_shift
      use caspt2_gradient, only: do_grad, do_csf, iRoot1, iRoot2
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"

#include "caspt2_grad.fh"
#include "pt2_guga.fh"
#include "chocaspt2.fh"
      DIMENSION DMAT(*),UEFF(nState,nState)
      Dimension VECROT(nState)


      IF (do_grad) THEN
        !! Print out some information for the first time only
    !       If (.not.IFMSCOUP.or.(IFMSCOUP.and.jState.eq.1))
    !  *      Call GradStart
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
        CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
C First, put in the reference density matrix.
        IDMOFF=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NORB(ISYM)
          DO II=1,NI
            IDM=IDMOFF+(II*(II+1))/2
            DMAT(IDM)=2.0D0
          END DO
          DO IT=1,NA
            ITABS=NAES(ISYM)+IT
            ITTOT=NI+IT
            DO IU=1,IT
              IUABS=NAES(ISYM)+IU
              IUTOT=NI+IU
              IDRF=(ITABS*(ITABS-1))/2+IUABS
              IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
              DMAT(IDM)=WORK(LDREF-1+IDRF)
            END DO
          END DO
           IDMOFF=IDMOFF+(NO*(NO+1))/2
        END DO
*       write(6,*)' DENS. Initial DMAT:'
*       WRITE(*,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
        CALL GETMEM('DPT','ALLO','REAL',LDPT,NDPT)
        CALL GETMEM('DSUM','ALLO','REAL',LDSUM,NDPT)
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDSUM),1)
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
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
        IF (IPRGLB.GE.verbose) THEN
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(6,'(a,2f10.2)')" Lambda  : CPU/WALL TIME=", cput,wallt
        END IF
C
C
C
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        !! Diagonal part
        CALL TRDNS2D(iVecX,iVecR,WORK(LDPT),NDPT,VECROT(JSTATE))
        CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*       write(6,*)' DPT after TRDNS2D.'
*       WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
        !! Off-diagonal part, if full-CASPT2
        IF (MAXIT.NE.0) THEN
          !! off-diagonal are ignored for CASPT2-D
          CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
          CALL TRDNS2O(iVecX,iVecR,WORK(LDPT),VECROT(JSTATE))
          CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
        END IF
*       write(6,*)' DPT after TRDNS2O.'
*       WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
        CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
        IF (IPRGLB.GE.verbose) THEN
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(6,'(a,2f10.2)')" TRDNS2DO: CPU/WALL TIME=", cput,wallt
        END IF
C
        !! D^PT in MO
        CALL GETMEM('DPT   ','ALLO','REAL',ipDPT   ,nDPTAO)
        !! D^PT(C) in MO
        CALL GETMEM('DPTC  ','ALLO','REAL',ipDPTC  ,nDPTAO)
        !! DPTAO1 (D^PT in AO, but not DPTA-01) couples with
        !! the CASSCF density (assume state-averaged) through ERIs.
        !! This density corresponds to the eigenvalue derivative.
        !! This is sometimes referred to as DPT2(AO) else where.
        CALL GETMEM('DPTAO ','ALLO','REAL',ipDPTAO ,nDPTAO)
        !! DPTAO2 couples with the inactive density.
        !! This density comes from derivative of the generalized
        !! Fock matrix (see for instance Eq. (24) in the 1990 paper).
        !! This is sometimes referred to as DPT2C(AO) else where.
        CALL GETMEM('DPTCAO','ALLO','REAL',ipDPTCAO,nDPTAO)
        !! DPTAO,DPTCAO,FPTAO,FPTCAO are in a block-squared form
        CALL GETMEM('FPT   ','ALLO','REAL',ipFPT   ,nDPTAO)
        CALL GETMEM('FPTC  ','ALLO','REAL',ipFPTC  ,nDPTAO)
        CALL GETMEM('FPTAO ','ALLO','REAL',ipFPTAO ,nDPTAO)
        CALL GETMEM('FPTCAO','ALLO','REAL',ipFPTCAO,nDPTAO)
        !! Transformation matrix
        CALL GETMEM('TRFMAT','ALLO','REAL',ipTrf   ,nBsqT)
        nch=0
        If (IfChol) nch=nvloc_chobatch(1)
        CALL GETMEM('WRK1  ','ALLO','REAL',ipWRK1  ,Max(nBasT**2,nch))
        CALL GETMEM('WRK2  ','ALLO','REAL',ipWRK2  ,Max(nBasT**2,nch))
        !! FIFA and FIMO (due to frozen orbitals)
C       CALL GETMEM('FIFA  ','ALLO','REAL',ipFIFA  ,nBsqT)
C       CALL GETMEM('FIMO  ','ALLO','REAL',ipFIMO  ,nBsqT)
        !! state-averaged density
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA ,nAshT*nAshT)
        !! Derivative of state-averaged density
        CALL GETMEM('RDMEIG','ALLO','REAL',ipRDMEIG,nAshT*nAshT)
C       write(6,*) "olag before"
C       call sqprt(Work(ipolag),nbast)
C
        Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPT),1)
        Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPTC),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipDPTAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipDPTCAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPT),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTC),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTCAO),1)
        If (.not.IfChol) Then
        Call DCopy_(nBsqT ,[0.0D+00],0,Work(ipFIFA),1)
        Call DCopy_(nBsqT ,[0.0D+00],0,Work(ipFIMO),1)
        End If
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipRDMSA),1)
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipRDMEIG),1)
C
        Call DCopy_(nCLag,[0.0D+00],0,Work(ipCLag),1)
        Call DCopy_(nOLag,[0.0D+00],0,Work(ipOLag),1)
C
        If (nFroT.ne.0) Then
          CALL GETMEM('DIA   ','ALLO','REAL',ipDIA ,nBsqT)
          CALL GETMEM('DI    ','ALLO','REAL',ipDI  ,nBsqT)
        End If
C
        If (do_csf) Then
          CALL GETMEM('DPTCanti','ALLO','REAL',ipDPTCanti,nDPTAO)
          Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPTCanti),1)
        End If
C
        !! Work(LDPT) -> Work(ipDPT2)
        !! Note that Work(ipDPT2) has the index of frozen orbitals.
        !! Note also that unrelaxed (w/o Z-vector) dipole moments with
        !! frozen orbitals must be wrong.
C       call dcopy_(ndpt,[0.0d+00],0,work(ldpt),1)
        If (nFroT.eq.0) Then
          Call DCopy_(nOsqT,Work(LDSUM),1,Work(ipDPT),1)
        Else
          Call OLagFro0(Work(LDSUM),Work(ipDPT))
        End If
C
        !! Construct the transformation matrix
        !! It seems that we have to transform quasi-canonical
        !! to CASSCF orbitals. The forward transformation has been
        !! done in ORBCTL.
        !!   C(PT2) = C(CAS)*X    ->    C(CAS) = C(PT2)*X^T
        !!   -> L(CAS) = X*L(PT2)*X^T
        !! inactive and virtual orbitals are not affected.
        Call DCopy_(nBsqT,[0.0D+0],0,Work(ipTrf),1)
        Call CnstTrf(Work(LTOrb),Work(ipTrf))
C       call sqprt(work(iptrf),nbast)
C
        !! Construct state-averaged density matrix
        Call DCopy_(nDRef,[0.0D+00],0,Work(ipWRK1),1)
        If (IFSADREF) Then
          Do iState = 1, nState
            Wgt  = 1.0D+00/nState
            Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
     *                  Work(ipWRK1),1)
          End Do
        Else
          Call DaXpY_(nDRef,1.0D+00,Work(LDMix+nDRef*(jState-1)),1,
     *                Work(ipWRK1),1)
        End If
        Call SQUARE(Work(ipWRK1),Work(ipRDMSA),1,nAshT,nAshT)
C       write(6,*) "state-averaged density matrix"
C       call sqprt(work(iprdmsa),nasht)
C
C       ----- Construct configuration Lagrangian -----
C
        !! For CI coefficient derivatives (CLag)
        !! Calculate the configuration Lagrangian
        !! This is done in the quasi-canonical basis
        CALL GETMEM('DEPSA ','ALLO','REAL',ipDEPSA,nAshT*nAshT)
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipDEPSA),1)
        !! Derivative of off-diagonal H0 of <Psi1|H0|Psi1>
        IF (MAXIT.NE.0) Call SIGDER(iVecX,iVecR,VECROT(jState))
        Call CLagX(1,Work(ipCLag),Work(ipDEPSA),VECROT)
C       call test3_dens(work(ipclag))
C       write(6,*) "original depsa"
C       call sqprt(work(ipdepsa),nasht)
C       write(6,*) "original depsa (sym)"
          do i = 1, nasht
          do j = 1, i-1
            val =(work(ipdepsa+i-1+nasht*(j-1))
     *           +work(ipdepsa+j-1+nasht*(i-1)))*0.5d+00
            work(ipdepsa+i-1+nasht*(j-1)) = val
            work(ipdepsa+j-1+nasht*(i-1)) = val
          end do
          end do
C       call sqprt(work(ipdepsa),nasht)
C
        If (NRAS1T+NRAS3T.NE.0) Then
          !! The density of the independent pairs (off-diagonal blocks)
          !! should be determined by solving Z-vector, so these blocks
          !! should be removed...?
C         write(6,*) "removing DEPSA of off-diagonal blocks"
C         write(6,*) "before"
C         call sqprt(work(ipdepsa),nasht)
            Do II = 1, nRAS1T
              Do JJ = nRAS1T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
            Do II = nRAS1T+1, nRAS1T+nRAS2T
              Do JJ = nRAS1T+nRAS2T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
C         write(6,*) "after"
C         call sqprt(work(ipdepsa),nasht)
          IF (IPRGLB.GE.debug)
     *      write(6,*) "depsa (sym) after removing off-diagonal blocks"
        Else
          IF (IPRGLB.GE.debug)
     *      write(6,*) "depsa (sym)"
        End If
        IF (IPRGLB.GE.verbose) call sqprt(work(ipdepsa),nasht)
C
        !! Configuration Lagrangian for MS-CASPT2
        !! This is the partial derivative of the transition reduced
        !! density matrices
        If (IFMSCOUP) Then
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          Call DerHEff(Work(ipCLag),VECROT)
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB.GE.verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(6,'(a,2f10.2)')" DerHEff : CPU/WALL TIME=", cput,wallt
            write(6,*)
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
        If (.not.INVAR) Then
          !! But, save the diagonal elements
          CALL GETMEM('DEPSA ','ALLO','REAL',ipDEPSAD,nAshT)
          Call DCopy_(nAshT,Work(ipDEPSA),nAshT+1,Work(ipDEPSAD),1)
          !! Clear
          Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipDEPSA),1)
        End If
C       write(6,*) "depsad"
C       call sqprt(work(ipdepsa),nasht)
C
        !! Transform the quasi-variational amplitude (T+\lambda/2?)
        !! in SR (iVecX) to C (iVecC2)
        !! Note that the contribution is multiplied by two
        !! somewhere else (maybe in olagns?)
        If (real_shift .ne. 0.0D+00 .or. imag_shift .ne. 0.0D+00
     &      .OR. IFMSCOUP) Then
          !! Have to weight the T-amplitude for MS-CASPT2
          IF (IFMSCOUP) THEN
            !! add lambda
            CALL PLCVEC(VECROT(jState),0.50d+00,IVECX,IVECR)
            CALL PTRTOC(1,IVECR,IVECC2)
            !! T-amplitude
            Do iStLag = 1, nStLag
              If (iStLag.eq.jState) Cycle
              Scal = VECROT(iStLag)
              If (ABS(Scal).LE.1.0D-12) Cycle
              Call MS_Res(2,jStLag,iStLag,Scal*0.5d+00)
            End Do
            If (do_csf) Then
              !! Prepare for something <\Phi_K^{(1)}|Ers|L>
              Call RHS_ZERO(7)
              ibk = ivecc2
              ivecc2 = 7
              Do iStLag = 1, nStLag
                If (iStLag.eq.jState) Cycle
                Scal = UEFF(iStLag,iRoot1)*UEFF(jStLag,iRoot2)
     *               - UEFF(jStLag,iRoot1)*UEFF(iStLag,iRoot2)
                Scal = Scal*0.5d+00
                If (ABS(Scal).LE.1.0D-12) Cycle
                Call MS_Res(2,jStLag,iStLag,Scal)
              End Do
              ivecc2 = ibk
            End If
          ELSE
            !! Add lambda to the T-amplitude
            CALL PLCVEC(0.5D+00,1.0D+00,IVECR,IVECX)
            CALL PTRTOC(1,IVECX,IVECC2)
          END IF
        End If
C
C         ipTrfL = ipTrf+nAshT*nBasT+nAshT
C         Call DGemm_('n','N',nAshT,nAshT,nAshT,
C    *                1.0D+00,Work(ipTrfL),nBasT,Work(ipDEPSA),nAshT,
C    *                0.0D+00,Work(ipdptcao),nAshT)
C         Call DGemm_('N','t',nAshT,nAshT,nAshT,
C    *                1.0D+00,Work(ipdptcao),nAshT,Work(ipTrfL),nBasT,
C    *                0.0D+00,Work(ipDEPSA),nAshT)
C
C       !! Just add DEPSA to DPT2
        Call AddDEPSA(Work(ipDPT),Work(ipDEPSA))
        !! Just transform the density in MO to AO
        CALL DPT2_Trf(Work(LDPT),Work(ipDPTAO),
     *                Work(LCMOPT2),Work(ipDEPSA),
     *                Work(LDSUM))
C       CALL GETMEM('DEPSA ','FREE','REAL',ipDEPSA,nAshT)
        !! Save the AO density
        !! ... write
C
C       ----- Construct orbital Lagrangian -----
C
        If (nFroT.ne.0) Then
          !! If frozen orbitals exist, we need to obtain
          !! electron-repulsion integrals with frozen orbitals to
          !! construct the orbital Lagrangian.
          If (.not.IfChol) Call TRAFRO(1)
C
          !! Get density matrix (Work(ipDIA)) and inactive density
          !! matrix (Work(ipDI)) to compute FIFA and FIMO.
          Call OLagFroD(Work(ipDIA),Work(ipDI),Work(ipRDMSA),
     *                  Work(ipTrf))
C         write(6,*) "density matrix"
C         call sqprt(work(ipdia),12)
C         call sqprt(work(ipdi),12)
        End If
C
        !! Construct orbital Lagrangian that comes from the derivative
        !! of ERIs. Also, do the Fock transformation of the DPT2 and
        !! DPT2C densities.
        If (IfChol) Then
          NumChoTot = 0
          Do iSym = 1, nSym
            NumChoTot = NumChoTot + NumCho_PT2(iSym)
          End Do
          Call GetMem('A_PT2 ','ALLO','REAL',ipA_PT2,NumChoTot**2)
          Call dcopy_(NumChoTot**2,[0.0D+00],0,Work(ipA_PT2),1)
        End If
        Do iSym = 1, nSym
          nOcc  = nIsh(iSym)+nAsh(iSym)
          If (.not.IfChol.or.iALGO.ne.1) Then
            lT2AO = nOcc*nOcc*nBasT*nBasT
            Call GetMem('T2AO','Allo','Real',ipT2AO,lT2AO)
            Call DCopy_(lT2AO,[0.0D+00],0,Work(ipT2AO),1)
          End If
C
          !! Orbital Lagrangian that comes from the derivative of ERIs.
          !! OLagNS computes only the particle orbitals.
C         write(6,*) "ialgo = ", ialgo
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          If (IfChol.and.iALGO.eq.1) Then
            CALL OLagNS_RI(iSym,Work(ipDPTC),Work(ipDPTCanti),
     *                     Work(ipA_PT2),NumChoTot)
          Else
            CALL OLagNS2(iSym,Work(ipDPTC),Work(ipT2AO))
          End If
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB.GE.verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(6,'(a,2f10.2)')" OLagNS  : CPU/WALL TIME=", cput,wallt
          END IF
C         write(6,*) "DPT2C"
C         call sqprt(work(ipdptc),nbast)
C
          !! MO -> AO transformations for DPT2 and DPT2C
          If ((.not.IfChol.or.iALGO.ne.1).or.nFroT.eq.0) Then
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPT),
     *                   Work(ipDPTAO),Work(ipWRK1))
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPTC),
     *                   Work(ipDPTCAO),Work(ipWRK1))
C           write(6,*) "dpt2"
C           call sqprt(work(ipdpt),nbast)
C           write(6,*) "dpt2ao"
C           call sqprt(work(ipdptao),nbast)
          End If
C
          !! Do some transformations relevant to avoiding (VV|VO)
          !! integrals. Orbital Lagrangian for the hole orbitals are
          !! computed. At the same time, F = G(D) transformations are
          !! also performed for D = DPT2 and DPT2C
          !! The way implemented (what?) is just a shit. I cannot find
          !! FIFA and FIMO for frozen orbitals, so I have to construct
          !! them. Here is the transformation of G(D^inact) and G(D).
          !! Work(ipFIFA) and Work(ipFIMO) computed in this subroutine
          !! is not yet correct. They are just two-electron after this
          !! subroutine.
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          CALL OLagVVVO(iSym,Work(ipDPTAO),Work(ipDPTCAO),
     *                  Work(ipFPTAO),Work(ipFPTCAO),Work(ipT2AO),
     *                  Work(ipDIA),Work(ipDI),Work(ipFIFA),
     *                  Work(ipFIMO),Work(ipA_PT2),NumChoTot)
        !   write(6,*) "olag after vvvo"
        !   call sqprt(work(ipolag),nbast)
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB.GE.verbose) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(6,'(a,2f10.2)')" OLagVVVO: CPU/WALL TIME=", cput,wallt
          END IF
C     write(6,*) "OLag"
C     do i = 1, 144
C       write(6,'(i3,f20.10)') i,work(ipolag+i-1)
C     end do
C      write(6,*) "fpt2ao"
C     call sqprt(work(ipfptao),12)
C     call abend
C
          !! AO -> MO transformations for FPT2AO and FPT2CAO
          If ((.not.IfChol.or.iALGO.ne.1).or.nFroT.eq.0) Then
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPT),
     *                   Work(ipFPTAO),Work(ipWRK1))
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPTC),
     *                   Work(ipFPTCAO),Work(ipWRK1))
          End If
C
          If (.not.IfChol.or.iALGO.ne.1) Then
            Call GetMem('T2AO','Free','Real',ipT2AO,lT2AO)
          End If
        End Do
        If (IfChol) Then
          Call GetMem('A_PT2 ','FREE','REAL',ipA_PT2,NumChoTot**2)
        End If
        !! Add DPTC to DSUM for the correct unrelaxed density
        !! Also, symmetrize DSUM
        Call AddDPTC(Work(ipDPTC),Work(LDSUM))
C
c
c
C       Call SQUARE(Work(LFIFA),Work(ipFIFA),1,12,12)
C       CALL DGEMM_('N','T',12,12,12,
C    *              2.0D+00,work(ipfifa),12,work(ipdpt),12,
C    *              1.0D+00,Work(ipOLAG),12)
C          call test2_dens(work(ipolag),work(ipdepsa))
C
c
c
C       write(6,*) "fptao after olagns"
C       call sqprt(Work(ipfptao),nbast)
C       write(6,*) "fptcao after olagns"
C       call sqprt(Work(ipfptcao),nbast)
C       write(6,*) "olag after olagns"
C       call sqprt(Work(ipolag),nbast)
C
        !! If frozen orbitals exist, frozen-inactive part of the
        !! unrelaxed PT2 density matrix is computed using the orbital
        !! Lagrangian. Additionally, Fock transformation is also
        !! required.
        If (nFroT.ne.0) Then
          !! Compute DPT2 density for frozen-inactive
C         write(6,*) "dpt before frozen"
C         call sqprt(work(ipdpt),nbast)
          if (.not.ifchol) then
            !! Construct FIFA and FIMO
            Call OLagFro3(Work(ipFIFA),Work(ipFIMO),Work(ipWRK1),
     *                    Work(ipWRK2))
          end if
          !! Add explicit FIMO and FIFA contributions. Implicit
          !! contributions are all symmetric in frozen + inactive
          !! orbitals, so they do not contribute to frozen density
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIFA),nBasT,Work(ipDPT),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIFA),nBasT,Work(ipDPT),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          !! Save DPT in order to subtract later
          Call DCopy_(nDPTAO,Work(ipDPT),1,Work(ipWRK1),1)
C
          !! Now, compute pseudo-density using orbital Lagrangian
          Call OLagFro1(Work(ipDPT),Work(ipOLag))
C
          !! Subtract the orbital Lagrangian added above.
          !! It is computed again in EigDer
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *               -1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *               -1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *               -1.0D+00,Work(ipFIFA),nBasT,Work(ipWRK1),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *               -1.0D+00,Work(ipFIFA),nBasT,Work(ipWRK1),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
C         write(6,*) "dpt after frozen"
C         call sqprt(work(ipdpt),nbast)
C
          !! Fock transformation for frozen-inactive density
          If (IfChol) Then
            iSym=1
            !! MO -> AO transformations for DPT2 and DPT2C
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPT),
     *                   Work(ipDPTAO),Work(ipWRK1))
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPTC),
     *                   Work(ipDPTCAO),Work(ipWRK1))
            !! For DF-CASPT2, Fock transformation of DPT2, DPT2C, DIA,
            !! DA is done here, but not OLagVVVO
            !! It seems that it is not possible to do this
            !! transformation in OLagVVVO, because the frozen-part of
            !! the DPT2 is obtained after OLagVVVO.
            Call OLagFro4(1,1,1,1,1,
     *                    Work(ipDPTAO),Work(ipDPTCAO),Work(ipFPTAO),
     *                    Work(ipFPTCAO),Work(ipWRK1))
            !! AO -> MO transformations for FPT2AO and FPT2CAO
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPT),
     *                   Work(ipFPTAO),Work(ipWRK1))
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPTC),
     *                   Work(ipFPTCAO),Work(ipWRK1))
          Else
C           write(6,*) "dpt"
C           call sqprt(work(ipdpt),nbast)
            Call OLagFro2(Work(ipDPT),Work(ipFPT),Work(ipWRK1),
     *                    Work(ipWRK2))
C         write(6,*) "fpt"
C           call sqprt(work(ipdpt),nbast)
          End If
C         write(6,*) "ipfifa"
C         call sqprt(work(ipfifa),12)
C         write(6,*) "ipfimo"
C         call sqprt(work(ipfimo),12)
      !   !! Construct FIFA and FIMO
      !   Call OLagFro3(Work(ipFIFA),Work(ipFIMO),Work(ipWRK1),
     *!                 Work(ipWRK2))
C
          CALL GETMEM('DIA   ','FREE','REAL',ipDIA ,nBsqT)
          CALL GETMEM('DI    ','FREE','REAL',ipDI  ,nBsqT)
        Else ! there are no frozen orbitals
          iSQ = 0
          iTR = 0
          Do iSym = 1, nSym
            nOrbI = nOrb(iSym)
            Call SQUARE(Work(LFIFA+iTR),Work(ipFIFA+iSQ),1,nOrbI,nOrbI)
            Call SQUARE(Work(LFIMO+iTR),Work(ipFIMO+iSQ),1,nOrbI,nOrbI)
            iSQ = iSQ + nOrbI*nOrbI
            iTR = iTR + nOrbI*(nOrbI+1)/2
          End Do
        End If
C         write(6,*) "ipfifa in dens"
C         call sqprt(work(ipfifa),nbast)
C         write(6,*) "ipfimo in dens"
C         call sqprt(work(ipfimo),nbast)
C        write(6,*) "FIFA in quasi-canonical"
C         call sqprt(work(ipfifa),12)
C        write(6,*) "FIFA in natural"
C         Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipTrf),nBasT,Work(ipFIFA),nBasT,
C    *                0.0D+00,Work(ipWRK1),nBasT)
C         Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *                0.0D+00,Work(ipWRK2),nBasT)
C         call sqprt(work(ipwrk2),12)
C
        !! Do some post-process for the contributions that comes from
        !! the above two densities.
C       CALL EigDer(Work(LDPT),Work(ipDPTC),Work(ipFPTAO),
C       write(6,*) "olag before eigder"
C       call sqprt(Work(ipolag),nbast)
C       write(6,*) "fpt2"
C       call sqprt(Work(ipfpt),nbast)
        CALL EigDer(Work(ipDPT),Work(ipDPTC),Work(ipFPTAO),
     *              Work(ipFPTCAO),Work(ipRDMEIG),Work(LCMOPT2),
     *              Work(ipTrf),Work(ipFPT),Work(ipFPTC),
     *              Work(ipFIFA),Work(ipFIMO),Work(ipRDMSA))
C          call test2_dens(work(ipolag),work(ipdepsa))
C       write(6,*) "olag after eigder"
C       call sqprt(Work(ipolag),nbast)
C       write(6,*) "Wlag after eigder"
C       call sqprt(work(ipwlag),nbast)
C       write(6,*) "rdmeig"
C       call sqprt(work(iprdmeig),nasht)
C       call abend
C
        !! Calculate the configuration Lagrangian again.
        !! The contribution comes from the derivative of eigenvalues.
        !! It seems that TRACI_RPT2 uses CI coefficients of RASSCF,
        !! so canonical -> natural transformation is required.
C       ipTrfL = ipTrf+nAshT*nBasT+nAshT
C       Call DGemm_('N','N',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipTrfL),nBasT,Work(ipRDMEIG),nAshT,
C    *              0.0D+00,Work(ipWRK1),nAshT)
C       Call DGemm_('N','T',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipWRK1),nAshT,Work(ipTrfL),nBasT,
C    *              0.0D+00,Work(ipRDMEIG),nAshT)
        If (.not.INVAR) Then !test
          CALL GETMEM('CLagT','ALLO','REAL',ipCLagT,nConf*nState)
          CALL GETMEM('EigT ','ALLO','REAL',ipEigT ,nAshT**2)
          Call DCopy_(nConf*nState,Work(ipCLag),1,Work(ipCLagT),1)
          Call DCopy_(nAshT**2,Work(ipRDMEIG),1,Work(ipEigT),1)
        End If
        !! Use canonical CSFs rather than natural CSFs in CLagEig
        ISAV = IDCIEX
        IDCIEX = IDTCEX
        !! Now, compute the configuration Lagrangian
        Call CLagEig(IFSSDM,Work(ipCLag),Work(ipRDMEIG))
C
        !! Now, here is the best place to compute the true off-diagonal
        !! active density for non-invariant CASPT2
        If (.not.INVAR) Then
          !! Add the density that comes from CI Lagrangian
          Call DEPSAOffC(Work(ipCLag),Work(ipDEPSA),Work(ipFIFA),
     *                   Work(ipFIMO),
     *                   Work(ipWRK1),Work(ipWRK2))
          !! Add the density that comes from orbital Lagrangian
          Call DEPSAOffO(Work(ipOLag),Work(ipDEPSA),Work(ipFIFA))
          !! Restore the diagonal elements
          Call DCopy_(nAshT,Work(ipDEPSAD),1,Work(ipDEPSA),nAshT+1)
          CALL GETMEM('DEPSAD','FREE','REAL',ipDEPSAD,nAshT)
          write(6,*) "DEPSA computed again"
          call sqprt(work(ipdepsa),nasht)
          If (NRAS1T+NRAS3T.NE.0) Then
            !! Remove the off-diagonal blocks for RASPT2
            Do II = 1, nRAS1T
              Do JJ = nRAS1T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
            Do II = nRAS1T+1, nRAS1T+nRAS2T
              Do JJ = nRAS1T+nRAS2T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
          End If
C         call dcopy_(nasht**2,[0.0d+00],0,work(ipdepsa),1)
C
          !! We have to do many things again...
          !! Just add DEPSA to DPT2
          Call AddDEPSA(Work(ipDPT),Work(ipDEPSA))
          !! Just transform the density in MO to AO
          CALL DPT2_Trf(Work(LDPT),Work(ipDPTAO),
     *                  Work(LCMOPT2),Work(ipDEPSA),
     *                  Work(LDSUM))
          !! Some transformations similar to EigDer
          Call EigDer2(Work(ipRDMEIG),Work(ipTrf),Work(ipFIFA),
     *                 Work(ipRDMSA),Work(ipDEPSA),
     *                 Work(ipWRK1),Work(ipWRK2))
C
          Call DCopy_(nConf*nState,Work(ipCLagT),1,Work(ipCLag),1) !test
         !test
          Call DaXpY_(nAshT**2,1.0D+00,Work(ipEigT),1,Work(ipRDMEIG),1)
          call DCopy_(nState*(nState-1)/2,[0.0D+00],0,Work(ipSLag),1)
          CALL GETMEM('CLagT','FREE','REAL',ipCLagT,nConf*nState)
          CALL GETMEM('EigT ','FREE','REAL',ipEigT ,nAshT**2)
          !! RDMEIG contributions
          !! Use canonical CSFs rather than natural CSFs
          !! Now, compute the configuration Lagrangian
          Call CLagEig(IFSSDM,Work(ipCLag),Work(ipRDMEIG))
          !! Now, compute the state Lagrangian and do some projections
          Call CLagFinal(Work(ipCLag),Work(ipSLag))
        End If
C
        !! Restore integrals without frozen orbitals, although not sure
        !! this operation is required.
        If (nFroT.ne.0.and..not.IfChol) Call TRAFRO(2)
C
        IDCIEX = ISAV
        !! Canonical -> natural transformation
        IF(ORBIN.EQ.'TRANSFOR') Then
          Do iState = 1, nState
            Call CLagX_TrfCI(Work(ipCLag+nConf*(iState-1)))
          End Do
        End If
        ! accumulate configuration Lagrangian only for MS,XMS,XDW,RMS,
        ! but not for SS-CASPT2
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          Call DaXpY_(nCLag,1.0D+00,Work(ipCLag),1,Work(ipCLagFull),1)
        end if
C       Call CLagFinal(Work(ipCLag),Work(ipSLag))
C
        !! Transformations of DPT2 in quasi-canonical to natural orbital
        !! basis and store the transformed density so that the MCLR
        !! module can use them.
        ! accumulate only if MS,XMS,XDW or RMS calculation
        ! call RecPrt('DPT2 before', '', work(ipDPT2), nBast, nBast)
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          Call DPT2_TrfStore(1.0D+00,Work(ipDPT),Work(ipDPT2),
     *                       Work(ipTrf),Work(ipWRK1))
          Call DPT2_TrfStore(2.0D+00,Work(ipDPTC),Work(ipDPT2C),
     *                       Work(ipTrf),Work(ipWRK1))
          If (do_csf) Then
            Call DPT2_TrfStore(1.0D+00,Work(ipDPTCanti),
     *      Work(ipDPT2Canti),Work(ipTrf),Work(ipWRK1))
          End If
        end if
        ! call RecPrt('DPT2 after', '', work(ipDPT2), nBast, nBast)
C       !! Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipTrf),nBasT,Work(LDPT),nBasT,
C    *              0.0D+00,Work(ipWRK1),nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *              0.0D+00,Work(ipWRK2),nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C         Call DaXpY_(nSQ,1.0D+00,Work(ipWRK2+iSQ),1,Work(ipDPT2+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C
C       !! Do the same for DPT2C Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipTrf),nBasT,Work(ipDPTC),nBasT,
C    *              0.0D+00,Work(ipWRK1),nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *              0.0D+00,Work(ipWRK2),nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C        Call DaXpY_(nSQ,2.0D+00,Work(ipWRK2+iSQ),1,Work(ipDPT2C+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C       call abend()
C       call sqprt(Work(ipRDMEIG),nAshT)
C
        !! square -> triangle so that the MCLR module can use the AO
        !! densities. Do this for DPT2AO and DPT2CAO (defined in
        !! caspt2_grad.f and caspt2_grad.h).
        ! accumulate only if MS,XMS,XDW or RMS calculation
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          iBasTr = 1
          iBasSq = 1
          Do iSym = 1, nSym
            nBasI = nBas(iSym)
            liBasTr = iBasTr
            liBasSq = iBasSq
            Do iBasI = 1, nBasI
              Do jBasI = 1, iBasI
                liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
                If (iBasI.eq.jBasI) Then
                  Work(ipDPT2AO +liBasTr-1) = Work(ipDPTAO +liBasSq-1)
                  Work(ipDPT2CAO+liBasTr-1) = Work(ipDPTCAO+liBasSq-1)
                Else
                  Work(ipDPT2AO +liBasTr-1)
     *            = Work(ipDPTAO +liBasSq-1)*2.0D+00
                  Work(ipDPT2CAO+liBasTr-1)
     *            = Work(ipDPTCAO+liBasSq-1)*2.0D+00
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
        If (IFSSDM.and.(jState.eq.iRlxRoot.or.nStLag.gt.1)) Then
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          If (.not.INVAR) Then
            write(6,*) "SS density matrix with IPEA not implemented"
            Call abend()
          End If
C
          !! Construct the SCF density
          !! We first need to construct the density averaged over all
          !! roots involved in SCF.
          Call DCopy_(nDRef,[0.0D+00],0,Work(ipWRK1),1)
          Call GetMem('LCI','ALLO','REAL',LCI,nConf)
          Wgt  = 1.0D+00/nState
          Do iState = 1, nState
C           Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
C    *                  Work(ipWRK1),1)
            Call LoadCI(WORK(LCI),iState)
            call POLY1(WORK(LCI))
            call GETDREF(WORK(ipWRK2))
            Call DaXpY_(nDRef,Wgt,Work(ipWRK2),1,Work(ipWRK1),1)
          End Do
          Call GetMem('LCI','FREE','REAL',LCI,nConf)
          !! Work(ipWRK2) is the SCF density (for nstate=nroots)
          Call SQUARE(Work(ipWRK1),Work(ipWRK2),1,nAshT,nAshT)
          Call DaXpY_(nAshT**2,-1.0D+00,Work(ipWRK2),1,Work(ipRDMSA),1)
          !! Construct the SS minus SA density matrix in Work(ipWRK1)
          Call OLagFroD(Work(ipWRK1),Work(ipWRK2),Work(ipRDMSA),
     *                  Work(ipTrf))
          !! Subtract the inactive part
          Call DaXpY_(nBasT**2,-1.0D+00,Work(ipWRK2),1,Work(ipWRK1),1)
          !! Here we should use ipDPTAO2??
          !! Save
          If (IfChol) Then
            Call CnstAB_SSDM(Work(ipDPTAO),Work(ipWRK1))
          Else
            !! Well, it is not working any more. I need to use
            !! Position='APPEND', but it is not possible if I need to
            !! use molcas_open or molcas_open_ext2
            write(6,*) "It is not possible to perform this calculation"
            write(6,*) "(non-state averaged density without"
            write(6,*) "density-fitting or Cholesky decomposition)"
            write(6,*) "Please use DF or CD"
            call abend()

          End If
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          IF (IPRGLB.GE.VERBOSE) THEN
            CPUT =CPTF10-CPTF0
            WALLT=TIOTF10-TIOTF0
            write(6,'(a,2f10.2)')" SSDM    : CPU/WALL TIME=", cput,wallt
          END IF
        End If
C       write(6,*) "pt2ao"
C       call sqprt(Work(ipDPTAO),12)
C       call prtril(Work(ipDPT2AO),12)
        CALL GETMEM('DEPSA ','FREE','REAL',ipDEPSA,nAshT)
C
        CALL GETMEM('DPT   ','FREE','REAL',ipDPT   ,nDPTAO)
        CALL GETMEM('DPTC  ','FREE','REAL',ipDPTC  ,nDPTAO)
        CALL GETMEM('DPTAO ','FREE','REAL',ipDPTAO ,nDPTAO)
        CALL GETMEM('DPTCAO','FREE','REAL',ipDPTCAO,nDPTAO)
        CALL GETMEM('FPT   ','FREE','REAL',ipFPT   ,nDPTAO)
        CALL GETMEM('FPTC  ','FREE','REAL',ipFPTC  ,nDPTAO)
        CALL GETMEM('FPTAO ','FREE','REAL',ipFPTAO ,nDPTAO)
        CALL GETMEM('FPTCAO','FREE','REAL',ipFPTCAO,nDPTAO)
C       CALL GETMEM('FIFA  ','FREE','REAL',ipFIFA  ,nBsqT)
C       CALL GETMEM('FIMO  ','FREE','REAL',ipFIMO  ,nBsqT)
        If (do_csf) Then
          CALL GETMEM('DPTCanti','FREE','REAL',ipDPTCanti,nDPTAO)
        End If
C
C       call test_dens(work(ipolag),work(ipclag),work(iptrf),
C    *                 work(ipwrk1),work(ipwrk2))
C
C       CALL GETMEM('WRK1  ','ALLO','REAL',ipWRK1,nBasT*nBasT)
C       CALL GETMEM('WRK2  ','ALLO','REAL',ipWRK2,nBasT*nBasT)
C       call dcopy_(nbast*nbast,[0.0d+00],0,Work(ipWRK1),1)
C       do i = 1, 5
C         Work(ipWRK1+i-1+nBasT*(i-1)) = 1.0D+00
C       end do
C       do i = 1, 5
C         ii = 5+i
C         do j = 1, 5
C           jj = 5+j
C           Work(ipWRK1+ii-1+nBasT*(jj-1)) = Work(LTORB+25+i-1+5*(j-1))
C         End Do
C       end do
C       do i = 11, 12
C         Work(ipWRK1+i-1+nBasT*(i-1)) = 1.0D+00
C       end do
C       write(6,*) "square transformation matrix"
C       call sqprt(work(ipwrk1),12)
C       write(6,*) "OLag before transformation"
C       call sqprt(work(ipolag),12)
C       Do iSym = 1, nSym
C       write(6,*) "olag before"
C       call sqprt(work(ipolag),nbast)
      if (.false.) then

        Call OLagFinal(Work(ipOLag),Work(ipTrf))

      else
C
        CALL GETMEM('WLAGL  ','ALLO','REAL',ipWLagL  ,nWLag)
        Call DCopy_(nWLag,[0.0d+00],0,Work(ipWLagL),1)
        Call DaXpY_(nWLag,0.5D+00,Work(ipOLag),1,Work(ipWLagL),1)
C       write(6,*) "Wlag square"
C       call sqprt(work(ipwlag),nbast)
C
        !! W(MO) -> W(AO) using the quasi-canonical orbitals
        !! No need to back transform to natural orbital basis
        Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(LCMOPT2),nBasT,Work(ipWLagL),nBasT,
     *              0.0D+00,Work(ipWRK1),nBasT)
        Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(ipWRK1),nBasT,Work(LCMOPT2),nBasT,
     *              0.0D+00,Work(ipWLagL),nBasT)
C
        !! square -> triangle for WLag(AO)
        Call DCopy_(nBasT*nBasT,Work(ipWLagL),1,Work(ipWRK1),1)
        iBasTr = 1
        iBasSq = 1
        Do iSym = 1, nSym
          nBasI = nBas(iSym)
          liBasTr = iBasTr
          liBasSq = iBasSq
          Do iBasI = 1, nBasI
            Do jBasI = 1, iBasI
              liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
              If (iBasI.eq.jBasI) Then
                Work(ipWLagL  +liBasTr-1) = Work(ipWRK1  +liBasSq-1)
              Else
              liBasSq2 = iBasSq + jBasI-1 + nBasI*(iBasI-1)
                Work(ipWLagL  +liBasTr-1)
     *            = Work(ipWRK1  +liBasSq-1)
     *            + Work(ipWRK1  +liBasSq2-1)
              End If
              liBasTr = liBasTr + 1
            End Do
          End Do
          iBasTr = iBasTr + nBasI*(nBasI+1)/2
          iBasSq = iBasSq + nBasI*nBasI
        End Do
        ! accumulate W Lagrangian only for MS,XMS,XDW,RMS,
        ! but not for SS-CASPT2
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
            Call DaXpY_(nBasSq,1.0D+00,Work(ipWLagL),1,Work(ipWLag),1)
        end if
        CALL GETMEM('WLAGL  ','FREE','REAL',ipWLagL  ,nWLag)
C
C
C
        !! Transform quasi-canonical -> natural MO basis
        !! orbital Lagrangian
C       write(6,*) "OLag before transformation"
C       call sqprt(work(ipolag),12)
C       write(6,*) "OLag after antisymmetrization"
C       call sqprt(work(ipolag),12)
        Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(ipTrf),nBasT,Work(ipOLag),nBasT,
     *              0.0D+00,Work(ipWRK1),nBasT)
        Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
     *              0.0D+00,Work(ipOLag),nBasT)
C       write(6,*) "OLag after transformation (pre)"
C       call sqprt(work(ipolag),12)
        !! sufficient only for active
        nBasI = nBas(1)
        Call DCopy_(nBasI**2,Work(ipOLag),1,Work(ipWRK1),1)
        Call DGeSub(Work(ipWRK1),nBas(1),'N',
     &              Work(ipWRK1),nBas(1),'T',
     &              Work(ipOLag),nBas(1),
     &              nBas(1),nBas(1))
        ! accumulate orbital Lagrangian only for MS,XMS,XDW,RMS,
        ! but not for SS-CASPT2
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          Call DaXpY_(nOLag,1.0D+00,Work(ipOLag),1,Work(ipOLagFull),1)
        end if
        ! call RecPrt('OLag    ','',work(ipOLag)    ,nBasT,nBasT)
        ! call RecPrt('OLagFull','',work(ipOLagFull),nBasT,nBasT)
      end if
C
        CALL GETMEM('TRFMAT','FREE','REAL',ipTRF   ,nBsqT)
        CALL GETMEM('WRK1  ','FREE','REAL',ipWRK1,nBasT*nBasT)
        CALL GETMEM('WRK2  ','FREE','REAL',ipWRK2,nBasT*nBasT)
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA ,nAshT*nAshT)
        CALL GETMEM('RDMEIG','FREE','REAL',ipRDMEIG,nAshT*nAshT)
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
        CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
C First, put in the reference density matrix.
        IDMOFF=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NORB(ISYM)
          DO II=1,NI
            IDM=IDMOFF+(II*(II+1))/2
            DMAT(IDM)=2.0D0
          END DO
          DO IT=1,NA
            ITABS=NAES(ISYM)+IT
            ITTOT=NI+IT
            DO IU=1,IT
              IUABS=NAES(ISYM)+IU
              IUTOT=NI+IU
              IDRF=(ITABS*(ITABS-1))/2+IUABS
              IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
              DMAT(IDM)=WORK(LDREF-1+IDRF)
            END DO
          END DO
           IDMOFF=IDMOFF+(NO*(NO+1))/2
        END DO
*       WRITE(6,*)' DENS. Initial DMAT:'
*       WRITE(6,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
        CALL GETMEM('DPT','ALLO','REAL',LDPT,NDPT)
        CALL GETMEM('DSUM','ALLO','REAL',LDSUM,NDPT)
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDSUM),1)

C The 1st order contribution to the density matrix
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
        CALL TRDNS1(IVEC,WORK(LDPT))
        CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*       WRITE(6,*)' DPT after TRDNS1.'
*       WRITE(6,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
        CALL TRDNS2D(IVEC,IVEC,WORK(LDPT),NDPT,1.0D+00)
        IF(IFDENS) THEN
C The exact density matrix evaluation:
          CALL TRDTMP(WORK(LDPT))
        ELSE
C The approximate density matrix evaluation:
          CALL TRDNS2A(IVEC,IVEC,WORK(LDPT))
        END IF
        CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*       WRITE(6,*)' DPT after TRDNS2D.'
*       WRITE(6,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
        CALL TRDNS2O(IVEC,IVEC,WORK(LDPT),1.0D+00)
        CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
        ! WRITE(6,*)' DPT after TRDNS2O.'
        ! WRITE(6,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      END IF
C
      CALL GETMEM('DPT','FREE','REAL',LDPT,NDPT)
      IDMOFF=0
      IDSOFF=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        DO IP=1,NO
          DO IQ=1,IP
            IDM=IDMOFF+(IP*(IP-1))/2+IQ
            IDSUM=IDSOFF+IP+NO*(IQ-1)
            DMAT(IDM)=DMAT(IDM)+WORK(LDSUM-1+IDSUM)
          END DO
        END DO
        IDMOFF=IDMOFF+(NO*(NO+1))/2
        IDSOFF=IDSOFF+NO**2
      END DO
      CALL GETMEM('DSUM','FREE','REAL',LDSUM,NDPT)
C Scale with 1/DENORM to normalize
      X=1.0D0/DENORM
      If (do_grad) X=1.0D+00
      CALL DSCAL_(NDMAT,X,DMAT,1)

CSVC: For true parallel calculations, replicate the DMAT array
C so that the slaves have the same density matrix as the master.
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        IF (.NOT.KING()) THEN
          CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
        END IF
        CALL GADSUM(DMAT,NDMAT)
      END IF
#endif

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      Subroutine CnstTrf(Trf0,Trf)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"

      Dimension Trf0(*),Trf(*)

      iSQ = 0
      iTOrb = 1 ! LTOrb
      ipTrfL = 0
C     write(6,*) "norbt = ",norbt
C     write(6,*) "nosqt = ", nosqt
C     write(6,*) "nbast = ", nbast
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
C       write(6,*) "nBasI",nBas(iSym)
C       write(6,*) "nOrbI",nOrb(iSym)
C       write(6,*) "nFroI",nFro(iSym)
C       write(6,*) "nIshI",nIsh(iSym)
C       write(6,*) "nAshI",nAsh(iSym)
C       write(6,*) "nSshI",nSsh(iSym)
C       write(6,*) "nDelI",nDel(iSym)
        nCor  = nFroI + nIshI
        nVir  = nSshI + nDelI
        ipTrfL = ipTrfL + iSQ
        !! frozen + inactive
C       Do iIsh = 1, nFroI + nIshI
C         Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = 1.0D+00
C       End Do
        !! frozen
        Do iIsh = 1, nFroI
          Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = 1.0D+00
        End Do
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
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
C
      DIMENSION DPT2(*),DEPSA(nAshT,nAshT)
C
    !   write(6,*) "DPT2MO before active-active contribution"
    !   call sqprt(dpt2,nbas(1)-ndel(1))
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2.gt.0) Then
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
     *             +DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)))*0.5D+00
              DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1)) = Val
              DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
    !   write(6,*) "DPT2MO after DEPSA"
    !   call sqprt(dpt2,nbas(1)-ndel(1))
C
      End Subroutine AddDEPSA
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AddDPTC(DPTC,DSUM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
C
      DIMENSION DPTC(*),DSUM(*)
C
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)!-nDel(iSym)
        If (nOrbI2.gt.0) Then
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
     *             +DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)))*0.5D+00
              DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1)) = Val
              DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
C
      End Subroutine AddDPTC
C
C-----------------------------------------------------------------------
C
      Subroutine TRAFRO(MODE)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      DIMENSION nFroTmp(8),nOshTmp(8),nOrbTmp(8)
C
      If (Mode.eq.1) Then
        Do jSym = 1, 8
          nFroTmp(jSym) = nFro(jSym)
          nOshTmp(jSym) = nOsh(jSym)
          nOrbTmp(jSym) = nOrb(jSym)
          nOsh(jSym) = nFro(jSym)+nIsh(jSym)+nAsh(jSym)
          nOrb(jSym) = nOsh(jSym)+nSsh(jSym)
          nFro(jSym) = 0
        End Do
      End If
C
      Call GetMem('LCMO','ALLO','REAL',LCMO,NCMO)
      Call DCopy_(NCMO,WORK(LCMOPT2),1,WORK(LCMO),1)
      if (IfChol) then
        call TRACHO3(WORK(LCMO))
      else
        call TRACTL(0)
      end if
      Call GetMem('LCMO','FREE','REAL',LCMO,NCMO)
C
      If (Mode.eq.1) Then
        Do jSym = 1, 8
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
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
C
      Dimension DPT2q(*),DPT2n(*),Trf(*),WRK(*)
C
      iMO = 1
      Do iSym = 1, nSym
        If (nOrb(iSym).GT.0) Then
          nOrbI = nBas(iSym)-nDel(iSym)
          !! Quasi-canonical -> natural transformation of DPT2
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,
     *                1.0D+00,Trf(iMO),nOrbI,DPT2q(iMO),nOrbI,
     *                0.0D+00,WRK,nOrbI)
          Call DGemm_('N','T',nOrbI,nOrbI,nOrbI,
     *                   Scal,WRK,nOrbI,Trf(iMO),nOrbI,
     *                1.0D+00,DPT2n(iMO),nOrbI)
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
C
      Return
C
      End Subroutine DPT2_TrfStore
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DPT2_Trf(DPT2,DPT2AO,CMO,DEPSA,DSUM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      DIMENSION DPT2(*),DPT2AO(*),CMO(*),DEPSA(nAshT,nAshT),DSUM(*)
C
      !! DPT2 transformation
      !! Just transform DPT2 (in MO, block-squared) to DPT2AO (in AO,
      !! block-squared). Also, for DPT2C which couples with the inactive
      !! density matrix.
      CALL GETMEM('WRK   ','ALLO','REAL',ipWRK,nBSQT)
C
      !! MO -> AO back transformation
      iCMO =1
      iAO = 1
      iMO = 1
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
        If (nORB(ISYM).GT.0) Then
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
     *             +DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*0.5D+00
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
            End Do
          End Do
          !! First, DPT2 -> DPT2AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
     *                 1.0D+00,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,
     *                 0.0D+00,Work(ipWRK),nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
     *                 1.0D+00,Work(ipWRK),nBasI,CMO(iCMO),nBasI,
     *                 0.0D+00,DPT2AO(iAO),nBasI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
        iMO  = iMO  + nBasI*nBasI
      End Do
C
      CALL GETMEM('WRK   ','FREE','REAL',ipWRK,nBSQT)
C
      END SUBROUTINE DPT2_Trf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer(DPT2,DPT2C,FPT2AO,FPT2CAO,RDMEIG,CMO,Trf,
     *                  FPT2,FPT2C,FIFA,FIMO,RDMSA)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      DIMENSION DPT2(*),DPT2C(*),FPT2AO(*),FPT2CAO(*),RDMEIG(*),CMO(*),
     *          Trf(*)
      DIMENSION FPT2(*),FPT2C(*),FIFA(*),FIMO(*),RDMSA(*)
C
C     write (6,*) "here is eigder"
      CALL GETMEM('WRK1 ','ALLO','REAL',ipWRK1 ,nBSQT)
      CALL GETMEM('WRK2 ','ALLO','REAL',ipWRK2 ,nBSQT)
      CALL GETMEM('FPT2 ','ALLO','REAL',ipFPT2 ,nBSQT)
      CALL GETMEM('FPT2C','ALLO','REAL',ipFPT2C,nBSQT)
C
      !! AO -> MO transformation
      iCMO =1
      iAO = 1
C     iMO = 1
C     write(6,*) "fpt2ao"
C     call sqprt(fpt2ao,nbasT)
C     write(6,*) "fpt2cao"
C     call sqprt(fpt2cao,nbasT)
      if (nfrot.ne.0) then
        Call DCopy_(nBsqT,FPT2,1,Work(ipFPT2),1)
        Call DCopy_(nBsqT,FPT2C,1,Work(ipFPT2C),1)
      else
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
C       iOFF = iWTMP + nBas(iSym)*nBas(iSym)
        If (nOrb(iSym).GT.0) Then
          nBasI = nBas(iSym)
          nOrbI = nOrb(iSym)
          !! First, FPT2(AO) -> FPT2(MO)
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                 1.0D+00,CMO(iCMO),nBasI,FPT2AO(iAO),nBasI,
     *                 0.0D+00,Work(ipWRK1),nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                 1.0D+00,Work(ipWRK1),nOrbI,CMO(iCMO),nBasI,
     *                 0.0D+00,Work(ipFPT2+iAO-1),nOrbI)
          !! Second, FPT2C(AO) -> FPT2C(MO)
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                 1.0D+00,CMO(iCMO),nBasI,FPT2CAO(iAO),nBasI,
     *                 0.0D+00,Work(ipWRK1),nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                 1.0D+00,Work(ipWRK1),nOrbI,CMO(iCMO),nBasI,
     *                 0.0D+00,Work(ipFPT2C+iAO-1),nOrbI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
C       iMO  = iMO  + nBasI*nBasI
      End Do
      end if
C     write(6,*) "fpt2"
C     call sqprt(work(ipfpt2),nbasT)
C     write(6,*) "fpt2c"
C     call sqprt(work(ipfpt2c),nbasT)
C
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2) ,1)
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2C),1)
C     write(6,*) "fpt2mo"
C     call sqprt(work(ipfpt2),nbasT)
C
C     write(6,*) "fpt2 in MO"
C     do isym = 1, nsym
C       nbasi = nbas(isym)
C       write(6,*) "for symmetry :", isym,nbasi
C       call sqprt(work(ipfpt2),nbasi)
C     end do
C
      !! construct Fock in MO
C
C     write(6,*) "ndref = ", ndref
C     write(6,*) "nstate = ", state
C     write(6,*) "ldref"
C     do i = 1, ndref
C       write(6,'(i3,f20.10)') i,Work(ldref+i-1)
C     end do
C     write(6,*) "ldmix"
C     do istate = 1, nstate
C     write(6,*) "istate = ", istate
C     do i = 1, ndref
C       write(6,'(i3,f20.10)') i,Work(ldmix+i-1+ndref*(istate-1))
C     end do
C     end do
C     write(6,*) "average"
C     do i = 1, ndref
C       val = 0.0d+00
C       do istate = 1, nstate
C         val = val + Work(ldmix+i-1+ndref*(istate-1))/dble(nstate)
C       end do
C       write(6,'(i3,f20.10)') i,val
C     end do
      iSQ = 1
C     write(6,*) "olag before"
C     call sqprt(work(ipolag),nbast)
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        ! nSshI = nSsh(iSym)
        ! nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        CALL DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        !! Active orbital contributions: (p,q) = (all,act)
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA,nAshI*nAshI)
        !  Construct the active density of the orbital energy
        !  Assume the state-averaged density (SS- and XMS-CASPT2)
C       nSeq = 0
C       Call DCopy_(nAshI*nAshI,[0.0D+00],0,Work(ipWRK1),1)
C       Do iState = 1, nState
C         Wgt  = Work(LDWgt+iState-1+nState*(iState-1))
C         Wgt  = 1.0D+00/nState
C         Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
C    *                Work(ipWRK1),1)
C       End Do
        !  RDM of CASSCF
        !  This is likely defined by a set of natural orbitals.
        !  Here, we have to transform to a set of quasi-canonical
        !  orbitals, so forward transformation is appropriate.
C       Call SQUARE(Work(ipWRK1),Work(ipRDMSA),1,nAshI,nAshI)
        Call DCopy_(nAshT**2,RDMSA,1,Work(ipRDMSA),1)
        !! nbast?
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *                      Work(ipRDMSA),nAshT,
     *              0.0D+00,Work(ipWRK1),nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipWRK1),nAshT,
     *                      Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipRDMSA),nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              1.0D+00,Work(ipFPT2+iSQ-1+nOrbI*nCor),nOrbI,
     *                      Work(ipRDMSA),nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA,nAshI*nAshI)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
C       If (nFroI.eq.0) Then
C         Call SQUARE(Work(LFIFA+iSQ-1),Work(ipWRK1),1,nOrbI,nOrbI)
C       Else
C         Call OLagFroSq(iSym,Work(LFIFA+iSQ-1),Work(ipWRK1))
C       End If
C       write(6,*) "fock in MO"
C       call sqprt(FIFA(iSQ),norbi)
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
!    *              2.0D+00,Work(ipWRK1),nOrbI,DPT2(iSQ),nOrbI,
     *              2.0D+00,FIFA(iSQ),nOrbI,DPT2(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
C
        !! explicit derivative of the effective Hamiltonian
        !! dfpq/da = d/da(C_{mu p} C_{nu q} f_{mu nu})
        !!         = f_{mu nu}^a + (C_{mu m} U_{mp} C_{nu q}
        !!                       + C_{mu p} C_{nu m} U_{mq}) * f_{mu nu}
        !!         = f_{mu nu}^a + U_{mp} f_{mq} + U_{mq} f_{pm}
        !! U_{pq}  = f_{pm} df_{qm} + f_{mp} df_{mq}
C       Call SQUARE(Work(LFIMO+iSQ-1),Work(ipWRK1),1,nOrbI,nOrbI)
C       write(6,*) "effective fock in MO"
C       call sqprt(FIMO(iSQ),norbi)
C       If (nFroT.eq.0) Then
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
!    *              1.0D+00,Work(ipWRK1),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
        CALL DGEMM_('T','N',nOrbI,nOrbI,nOrbI,
!    *              1.0D+00,Work(ipWRK1),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
C       End If
        !! Implicit derivative of inactive orbitals (DPT2C)
        Call DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2C+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
C     write(6,*) "olag in eigder"
C     call sqprt(work(ipolag),nbasT)
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
     *        = RDMEIG(iSQA+iTUA) + Work(ipFPT2+iSq-1+iTU)
C           write(6,'(2i3,f20.10)') it,iu,Work(ipFPT2+iSq-1+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C     write(6,*) "rdmeig"
C     call sqprt(rdmeig,5)
C
      CALL GETMEM('WRK1 ','FREE','REAL',ipWRK1 ,nBSQT)
      CALL GETMEM('WRK2 ','FREE','REAL',ipWRK2 ,nBSQT)
      CALL GETMEM('FPT2 ','FREE','REAL',ipFPT2 ,nBSQT)
      CALL GETMEM('FPT2C','FREE','REAL',ipFPT2C,nBSQT)
C
      END SUBROUTINE EigDer
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer2(RDMEIG,Trf,FIFA,RDMSA,DEPSA,WRK1,WRK2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      DIMENSION RDMEIG(*),Trf(*)
      DIMENSION FIFA(*),RDMSA(*),DEPSA(*),WRK1(*),WRK2(*)
C
      CALL GETMEM('FPT2 ','ALLO','REAL',ipFPT2 ,nBSQT)
C
      !! Compute G(D), where D=DEPSA
      Call DEPSATrf(DEPSA,Work(ipFPT2),WRK1,WRK2)
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2),1)
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
        CALL DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        !! Active orbital contributions: (p,q) = (all,act)
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA,nAshI*nAshI)
        Call DCopy_(nAshT**2,RDMSA,1,Work(ipRDMSA),1)
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *                      Work(ipRDMSA),nAshT,
     *              0.0D+00,WRK1,nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,WRK1,nAshT,
     *                      Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipRDMSA),nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              1.0D+00,Work(ipFPT2+iSQ-1+nOrbI*nCor),nOrbI,
     *                      Work(ipRDMSA),nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA,nAshI*nAshI)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
        CALL DGEMM_('N','T',nOrbI,nAshI,nAshI,
     *              2.0D+00,FIFA(1+nOrbI*nCor),nOrbI,DEPSA,nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
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
            RDMEIG(iSQA+iTUA) = Work(ipFPT2+iSq-1+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C
      CALL GETMEM('FPT2 ','FREE','REAL',ipFPT2 ,nBSQT)
C
      END SUBROUTINE EigDer2
C
C-----------------------------------------------------------------------
C
      Subroutine DEPSATrf(DEPSA,FPT2,WRK1,WRK2)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      Dimension DEPSA(nAshT,nAshT),FPT2(*),WRK1(*),WRk2(*)
C
      Call DCopy_(nBasT**2,[0.0D+00],0,FPT2,1)
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
        Call GetMem('DAO ','ALLO','REAL',ipDAO ,nBsqT)
        Call GetMem('DMO ','ALLO','REAL',ipDMO ,nBsqT)
        Call GetMem('WRK1','ALLO','REAL',ipWRK1,nBsqT)
        Call GetMem('WRK2','ALLO','REAL',ipWRK2,nBsqT)
        !! First, MO-> AO transformation of DEPSA
        Do iSym = 1, nSym
          Call DCopy_(nBsqT,[0.0D+00],0,Work(ipDMO),1)
          nCorI = nFro(iSym)+nIsh(iSym)
          nBasI = nBas(iSym)
          Do iAsh = 1, nAsh(iSym)
            Do jAsh = 1, nAsh(iSym)
              Work(ipDMO+nCorI+iAsh-1+nBasI*(nCorI+jAsh-1))
     *          = DEPSA(iAsh,jAsh)
            End Do
          End Do
          Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDMO),
     *                 Work(ipDAO),Work(ipWRK1))
        End Do
        !! Compute G(D)
        Call DCopy_(nBsqT,[0.0D+00],0,Work(ipWRK1),1)
        Call DCopy_(nBsqT,[0.0D+00],0,Work(ipDMO),1)
        !! it's very inefficient
        Call OLagFro4(1,1,1,1,1,
     *                Work(ipDAO),Work(ipWRK1),Work(ipDMO),
     *                Work(ipWRK1),Work(ipWRK2))
        !! G(D) in AO -> G(D) in MO
        Do iSym = 1, nSym
          Call OLagTrf(2,iSym,Work(LCMOPT2),FPT2,
     *                 Work(ipDMO),Work(ipWRK1))
        End Do
        Call GetMem('DAO ','FREE','REAL',ipDAO ,nBsqT)
        Call GetMem('DMO ','FREE','REAL',ipDMO ,nBsqT)
        Call GetMem('WRK1','FREE','REAL',ipWRK1,nBsqT)
        Call GetMem('WRK2','FREE','REAL',ipWRK2,nBsqT)
      Else
        nCorI = nFro(iSym)+nIsh(iSym)
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            Call DaXpY_(nBasT**2,DEPSA(iAshI,jAshI),WRK1,1,FPT2,1)
C
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            Call DaXpY_(nBasT**2,-0.5D+00*DEPSA(iAshI,jAshI),
     *                  WRK1,1,FPT2,1)
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
C*MODULE MTHLIB  *DECK PRTRIL
      SUBROUTINE PRTRIL(D,N)
C
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION D(*)
C
      MAX = 5
      MM1 = MAX - 1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(6,9008)
         WRITE(6,9028) (I,I=I0,IL)
         WRITE(6,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            WRITE(6,9048) I,'        ',(D(J),J=J0,JL)
  100    CONTINUE
  120 CONTINUE
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,A8,10F11.6)
      END
C
C-----------------------------------------------------------------------
C
      Subroutine CnstAB_SSDM(DPT2AO,SSDM)
C
      use ChoVec_io
      use ChoSwp, only: InfVec
      use ChoArr, only: nDimRS
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
#include "warnings.h"
#include "chocaspt2.fh"
#include "choglob.fh"
C
      Dimension DPT2AO(*),SSDM(*)
      Character*4096 RealName
      Integer iSkip(8),ipWRK(8)
      integer nnbstr(8,3)
      Logical is_error
C
      ! INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      call getritrfinfo(nnbstr,maxvec,n2)
      iSym = 1 !! iSym0
C
      NumChoTot = 0
      Do jSym = 1, nSym
        NumChoTot = NumChoTot + NumCho_PT2(jSym)
      End Do
      NumCho=NumChoTot
      Do jSym = 1, nSym
        iSkip(jSym) = 1
      End Do
C
      nBasI  = nBas(iSym)
C
      Call GetMem('A_PT2 ','ALLO','REAL',ipA_PT2,NumChoTot**2)
      !! Read A_PT2
    !   Call PrgmTranslate('CMOPT2',RealName,lRealName)
    !   Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
    !  &                      'DIRECT','UNFORMATTED',
    !  &                      iost,.FALSE.,
    !  &                        1,'OLD',is_error)
    !   Read (LuCMOPT2) Work(ipA_PT2:ipA_PT2+NumChoTot**2-1)

      ! Read A_PT2 from LUAPT2
      id = 0
      call ddafile(LUAPT2, 2, work(ipA_PT2), numChoTot**2, id)

      !! Open B_PT2
      Call PrgmTranslate('GAMMA',RealName,lRealName)
      Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                      'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      nBas(iSym)**2*8,'OLD',is_error)
C
      CALL GETMEM('CHSPC','ALLO','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTVEC','ALLO','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','ALLO','REAL',ipWRK(iSym),nBasT*nBasT)
      !! V(P) = (mu nu|P)*D_{mu nu}
      CALL GETMEM('V1   ','ALLO','REAL',ipV1,NumCho)
      CALL GETMEM('V2   ','ALLO','REAL',ipV2,NumCho)
      !! B_SSDM(mu,nu,P) = D_{mu rho}*D_{nu sigma}*(rho sigma|P)
      Call GetMem('B_SSDM','ALLO','REAL',ipB_SSDM,NCHSPC)
C
      !! Prepare density matrix
      !! subtract the state-averaged density matrix
C
      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym).EQ.0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,iSym)
      JRED2=InfVec(NumCho_PT2(iSym),2,iSym)

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED.EQ.0) Cycle

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
          CALL CHO_VECRD(WORK(IP_CHSPC),NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)
          IF(NUMV.ne.JNUM) THEN
            write(6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(6,*)' read JNUM vectors. Instead it returned NUMV'
            write(6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC.NE.JRED) THEN
            write(6,*)' Rats! It was assumed that the Cholesky vectors'
            write(6,*)' in HALFTRNSF all belonged to a given reduced'
            write(6,*)' set, but they don''t!'
            write(6,*)' JRED, JREDC:',JRED,JREDC
            write(6,*)' Back to the drawing board?'
            write(6,*)' Let the program continue and see what happens.'
          END IF
C
          ipVecL = ip_CHSPC
          Do iVec = 1, NUMV
C
            !! reduced form -> squared AO vector (mu nu|iVec)
            ! If (l_NDIMRS.LT.1) Then
            If (size(nDimRS).lt.1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              ! lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
              lscr  = nDimRS(iSym,JREDL)
            End If
            Call DCopy_(nBasI**2,[0.0D+00],0,Work(ipWRK(iSym)),1)
            Call Cho_ReOrdr(irc,Work(ipVecL),lscr,1,
     *                      1,1,1,iSym,JREDC,2,ipWRK,
     *                      iSkip)
            ipVecL = ipVecL + lscr
C
            Work(ipV1+JV1-1+iVec-1) = DDot_(nBasI**2,DPT2AO,1,
     *                                      Work(ipWRK(iSym)),1)
            Work(ipV2+JV1-1+iVec-1) = DDot_(nBasI**2,SSDM  ,1,
     *                                      Work(ipWRK(iSym)),1)
C
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  1.0D+00,DPT2AO,nBasI,Work(ipWRK(iSym)),nBasI,
     *                  0.0D+00,Work(ipHTVec),nBasI)
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  1.0D+00,Work(ipHTVec),nBasI,SSDM,nBasI,
     *                  0.0D+00,Work(ipB_SSDM+nBasT**2*(iVec-1)),nBasI)
          End Do
          NUMVI = NUMV
C
          KV1=JSTART
          JBATCH_TOT=NBTCHES(iSym)
          DO JBATCH=1,NBATCH
            JBATCH_TOT=JBATCH_TOT+1

            KNUM=NVLOC_CHOBATCH(JBATCH_TOT)
            KV2=KV1+KNUM-1

            JREDC=JRED
            CALL CHO_VECRD(WORK(IP_CHSPC),NCHSPC,KV1,KV2,iSym,
     &                              NUMV,JREDC,MUSED)
           Call R2FIP(Work(ip_CHSPC),Work(ipWRK(iSym)),ipWRK(iSym),NUMV,
     *                size(nDimRS),infVec,nDimRS,
     *                nBasT,nSym,iSym,iSkip,irc,JREDC)
C
            !! Exchange part of A_PT2
            NUMVJ = NUMV
            CALL DGEMM_('T','N',NUMVI,NUMVJ,nBasT**2,
     *                 -1.0D+00,Work(ipB_SSDM),nBasT**2,
     *                          Work(ip_CHSPC),nBasT**2,
     *                1.0D+00,Work(ipA_PT2+JV1-1+NumCho*(KV1-1)),NumCho)
            KV1=KV1+KNUM
          END DO
C
          !! Read, add, and save the B_PT2 contribution
          Do iVec = 1, NUMVI
            Read  (Unit=LuGAMMA,Rec=JV1+iVec-1)
     *        Work(ipWRK(iSym):ipWRK(iSym)+nBasT**2-1)
            !! The contributions are doubled,
            !! because halved in PGet1_RI3?
            !! Coulomb
            Call DaXpY_(nBasT**2,Work(ipV2+JV1+iVec-2),
     *                  DPT2AO,1,Work(ipWRK(iSym)),1)
            Call DaXpY_(nBasT**2,Work(ipV1+JV1+iVec-2),
     *                  SSDM  ,1,Work(ipWRK(iSym)),1)
            !! Exchange
            Call DaXpY_(nBasT**2,-1.0D+00,
     *                  Work(ipB_SSDM+nBasT**2*(iVec-1)),1,
     *                  Work(ipWRK(iSym)),1)
            Write (Unit=LuGAMMA,Rec=JV1+iVec-1)
     *        Work(ipWRK(iSym):ipWRK(iSym)+nBasT**2-1)
          End Do
          JV1=JV1+JNUM
        End Do
      End Do
C
      Close (LuGamma)
C
      !! Coulomb for A_PT2
      !! Consider using DGER?
      Call DGEMM_('N','T',NumCho,NumCho,1,
     *            2.0D+00,Work(ipV1),NumCho,Work(ipV2),NumCho,
     *            1.0D+00,Work(ipA_PT2),NumCho)
C
      !! Write A_PT2
      ! REWIND LuCMOPT2
      ! Write (LuCMOPT2) Work(ipA_PT2:ipA_PT2+NumChoTot**2-1)
      ! Close (LuCMOPT2)

      ! write to A_PT2 in LUAPT2
      id = 0
      call ddafile(LUAPT2, 1, Work(ipA_PT2), NumChoTot**2, id)

      Call GetMem('A_PT2 ','FREE','REAL',ipA_PT2,NumChoTot**2)
C
      CALL GETMEM('CHSPC','FREE','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTVEC','FREE','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','FREE','REAL',ipWRK(iSym),nBasT*nBasT)
      CALL GETMEM('V1   ','FREE','REAL',ipV1,NumCho)
      CALL GETMEM('V2   ','FREE','REAL',ipV2,NumCho)
      Call GetMem('B_SSDM','FREE','REAL',ipB_SSDM,NCHSPC)
C     call abend
C
      End Subroutine CnstAB_SSDM
