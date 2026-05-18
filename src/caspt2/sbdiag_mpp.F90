!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

! SVC2010: The following subroutine computes the transformation
! matrices using global arrays rather than replicate arrays.  Currently,
! the maximum memory scales approximately as 3*(NAS**2), which is twice
! what is needed in the replicate routines.  It can be reduced to
! 2*(NAS**2), if the matrix multiplications are rewritten to occur in
! batch mode.  However, unlike in the replicate routine, this amount is
! divided over processors.
#ifdef _MOLCAS_MPP_
      SUBROUTINE SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
#ifdef _SCALAPACK_
      use scalapack_mod, only: GA_PDSYEVX_
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: INSANE
      USE Para_Info, ONLY: King
      use caspt2_global, only: LUSBT
      use caspt2_global, only: do_grad, do_lindep, nStpGrd, LUSTD,      &
     &                         idBoriMat
      use EQSOLV, only: IDSMAT, IDTMAT, IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nASup, nISup, Cases, IfDOrtho, ThrShn,   &
     &                         ThrShs, nInDep, BMATRIX, BTRANS, BSPECT
      use Constants, only: Zero, One
      use definitions, only: iwp, wp, u6
      IMPLICIT None

      integer(kind=iwp), intent(in):: iSym, iCase
      real(kind=wp), intent(out):: CondNr, CPU
!-SVC20100902: global arrays header files
#include "global.fh"
#include "mafdecls.fh"
#ifndef _SCALAPACK_
      real(kind=wp) WGRONK(2)
      real(kind=wp), allocatable:: VEC(:), SCRATCH(:)
      integer(kind=iwp) Info, NSCRATCH
#endif
      LOGICAL(kind=iwp) bSTAT
      CHARACTER(LEN=2) cSYM,cCASE
      real(kind=wp), allocatable:: COL(:), TMP(:), SD(:), SCA(:),       &
     &                             EIG(:), COND(:), TRANS(:), BD(:)
      integer(kind=iwp) NAS, NIS, NCOEF, lg_S, NCOL, NTMP, IOFF, J, IDS,&
     &                  MyRank, iLo, iHi, jLo, jHi, ISTA, IEND, MS, LDS,&
     &                  I, lg_V, NIN, mV, LDV, lg_T, IDT, lg_B, mB, lDB,&
     &                  IDB, IDB2, lg_X, lg_ST
      real(kind=wp) FP, SDiag, SZMIN, SZMAX, SZ, dTrans
      real(kind=wp) CPU1, CPUE, TIO, TIOE, CPU2
      real(kind=wp), External:: PSBMAT_FPRINT, DNRM2_

! On entry, the DRA metafiles contain the matrices S and B for cases A
! (iCASE=1) en C (iCASE=4).  These symmetric matrices are stored on disk
! in full square format.  The rectangular matrices T are computed, such
! that Sum_(J) [B(I,J)*T(J,MU)] = Sum_(J) [S(I,J)*T(J,MU)*BD(MU)], where
! I,J are in the range (1,NASUP(ISYM,ICASE)) and MU is in the range
! (1,NINDEP(ISYM,ICASE)).  NINDEP is the numerically effective rank of
! S, and the columns of T are orthonormal: Sum_(I,J)
! [S(I,J)*T(J,MU)*T(J,NU)] = Kron(MU,NU).  B is destroyed, and is
! overwritten by BD(MU) and T(I,MU), which is stored in a DRA metafile.

! Initialize the DRA I/O subsystem with default values.

      IF (iCASE.NE.1.AND.iCASE.NE.4) THEN
        WRITE(u6,*) 'Invalid CASE number used for global SBDIAG, Abort'
        CALL AbEnd()
      END IF

      write(unit=cCase, fmt='(I2.2)') iCase
      write(unit=cSYM, fmt='(I2.2)') iSYM
! Start a long loop over irreps:

      CPU=Zero
      CONDNR=Zero
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS
      IF(NCOEF.EQ.0) RETURN

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(u6,'("DEBUG> ",A12,A5,I2,A2,A6,A2,A5,I1)')                &
     &  'SBDIAG_MPP: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM
      END IF

! Allocate memory for the S matrix and read it from disk:
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)
      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_S,NAS)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ', FP
      END IF

! The S matrices are needed later on by non-global routines.  Take the
! oportunity to save them to LUSBT here.  FIXME: Should be removed once
! full parallelization of use of S matrices is achieved.
      IF (KING()) THEN
        NCOL=NAS
        NTMP=(NAS*(NAS+1))/2
        CALL mma_allocate(COL,NCOL,Label='COL')
        CALL mma_allocate(TMP,NTMP,Label='TMP')
        iOFF=0
        DO J=1,NAS
          call GA_Get (lg_S, 1, J, J, J, COL, NAS)
          CALL DCOPY_(J,COL,1,TMP(1+iOFF),1)
          iOFF=iOFF+J
        END DO
        CALL mma_deallocate(COL)
        IDS=IDSMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,TMP,NTMP,IDS)
        CALL mma_deallocate(TMP)
      END IF

! Save the diagonal elements from the S matrix for easy access later on.
! FIXME: nicer way to do this?
      CALL mma_allocate(SD,NAS,Label='SD')
      SD(:)=Zero
      myRank = GA_NodeID()
      call GA_Distribution (lg_S, myRank, iLo, iHi, jLo, jHi)
      ISTA=MAX(ILO,JLO)
      IEND=MIN(IHI,JHI)
      IF (ISTA.NE.0) THEN
        call GA_Access (lg_S, iLo, iHi, jLo, jHi, mS, LDS)
        DO I=ISTA,IEND
          SD(I)=DBL_MB(mS+I-ILO+LDS*(I-JLO))
        END DO
        call GA_Release (lg_S, iLo, iHi, jLo, jHi)
      END IF
      CALL GADGOP (SD,NAS,'+')

! Calculate the scaling factors and store them in array SCA.
      CALL mma_allocate(SCA,NAS,Label='SCA')
      DO I=1,NAS
        SDiag=SD(I)
        IF (IFDORTHO) THEN
          SCA(I)=One
        ELSE
          IF(SDiag.GT.THRSHN) THEN
            SCA(I)=(One+DBLE(I)*3.0E-6_wp)/SQRT(SDiag)
          ELSE
            SCA(I)=Zero
          END IF
        END IF
      END DO

! Scale the elements S(I,J) with the factor SCA(I)*SCA(J).
      myRank = GA_NodeID()
      call GA_Distribution (lg_S, myRank, iLo, iHi, jLo, jHi)
      IF (iLo.NE.0) THEN
        call GA_Access (lg_S, iLo, iHi, jLo, jHi, mS, LDS)
        call S_SCALE (NAS,SCA,DBL_MB(mS),iLo,iHi,jLo,jHi,LDS)
        call GA_Release_Update (lg_S, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()
      CALL TIMING(CPU1,CPUE,TIO,TIOE)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_S,NAS)
        WRITE(u6,'("DEBUG> ",A,ES21.14)')'SMAT NORM AFTER SCALING: ', FP
      END IF

      CALL mma_allocate(EIG,NAS,Label='EIG')
      EIG(:)=Zero
! Diagonalize the global array S.  Some (old) reports about parallel
! performance recommend PDSYEVX or PDSYEVR as fastest methods if
! eigenvectors are needed (FIXME: should time this).  For the linear
! dependence removal, split eigenvectors in horizontal stripes so that
! each processor has a row window of all column vectors
#ifdef _SCALAPACK_
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      CALL GA_PDSYEVX_ (lg_S, lg_V, EIG, 0)
      bSTAT = GA_Destroy (lg_S)
#else
! here for the non-ScaLAPACK version: copy matrix to master process,
! diagonalize using the serial DSYEV routine, and copy the resulting
! eigenvectors back to a global array.  Then distribute the eigenvalues.
      IF (myRank.EQ.0) THEN
        CALL mma_allocate(VEC,NAS**2,Label='VEC')
        CALL GA_Get (lg_S, 1, NAS, 1, NAS, VEC, NAS)
      END IF
      bSTAT = GA_Destroy (lg_S)
      IF (myRank.EQ.0) THEN
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        CALL mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,                            &
     &              SCRATCH,NSCRATCH,INFO)
        CALL mma_deallocate(SCRATCH)
      END IF
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      IF (myRank.EQ.0) THEN
        CALL GA_Put (lg_V, 1, NAS, 1, NAS, VEC, NAS)
        CALL mma_deallocate(VEC)
      END IF
      CALL GADGOP(EIG,NAS,'+')
#endif

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS,EIG,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'SMAT EIGENVALUE NORM: ', FP
      END IF

      NIN=0
      DO J=1,NAS
        IF(EIG(J).GE.THRSHS) NIN=NIN+1
      END DO
      NINDEP(ISYM,ICASE)=NIN
      IF (NIN.EQ.0) THEN
        CALL mma_deallocate(SCA)
        CALL mma_deallocate(EIG)
        CALL mma_deallocate(SD)
        bSTAT = GA_Destroy (lg_V)
        RETURN
      END IF

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1

      CALL mma_allocate(COND,NIN,Label='COND')
      COND(:)=Zero
! Form orthonormal transformation vectors by scaling the eigenvectors.
      call GA_Sync()
      myRank = GA_NodeID()
      call GA_Distribution (lg_V, myRank, iLo, iHi, jLo, jHi)
      IF (iLo.NE.0) THEN
        call GA_Access (lg_V, iLo, iHi, jLo, jHi, mV, LDV)
        IF ((jHi-jLo+1).NE.NAS) THEN
          WRITE(u6,*) 'SBDIAG_MPP: error in striping of lg_V, ABORT'
          CALL ABEND()
        END IF
        call V_SCALE (EIG,SCA(iLo),DBL_MB(mV),                          &
     &              iHi-iLo+1,jHi-jLo+1,LDV,NIN,COND)
        call GA_Release_Update (lg_V, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()

      CALL mma_deallocate(EIG)
      CALL mma_deallocate(SCA)

! The condition number, after scaling, disregarding linear dep.
! FIXME: adapt to local subroutine for global array lg_V
      IF(NIN.GE.2) THEN
        CALL GADGOP (COND,NIN,'+')
        SZMIN=1.0E99_wp
        SZMAX=Zero
        DO I=1,NIN
          SZ=COND(I)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=SZMAX/SZMIN
      END IF
      CALL mma_deallocate(COND)

! Copy the NIN non-linear dependent eigenvectors to the transformation
! matrix T(NAS,NIN).
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
      call GA_Copy_Patch ('N', lg_V, 1, NAS, 1, NIN,                    &
     &                         lg_T, 1, NAS, 1, NIN)
      bStat = GA_Destroy (lg_V)

      IF(BMATRIX.EQ.'NO      ') THEN
! In some calculations, we do not use B matrices.
! Write the T matrix to disk and exit.  FIXME: This
! should be removed when the transformation matrices are stored as disk
! resident arrays only.
        IF (KING()) THEN
          CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
          call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,TRANS,1)
            WRITE(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
          CALL mma_deallocate(TRANS)
        END IF
        bStat = GA_Destroy (lg_T)
        RETURN
      ELSE IF(BTRANS.NE.'YES') THEN
! In other calculations, the B matrix is used but not transformed.  We
! may need the diagonal active energies, i.e. the diagonal values of B
! divided by the diagonal values of S. These are placed where the
! eigenvalues would go in ordinary CASPT2.
        CALL PSBMAT_GETMEM ('B',lg_B,NAS)
        CALL PSBMAT_READ ('B',iCase,iSym,lg_B,NAS)
        CALL mma_allocate(BD,NAS,Label='BD')
        BD(:)=Zero
        myRank = GA_NodeID()
        call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
        ISTA=MAX(ILO,JLO)
        IEND=MIN(IHI,JHI)
        IF (ISTA.NE.0) THEN
          call GA_Access (lg_B, iLo, iHi, jLo, jHi, mB, LDB)
          DO I=ISTA,IEND
            BD(I)=DBL_MB(mB+I-ILO+LDB*(I-JLO))
          END DO
          call GA_Release (lg_B, iLo, iHi, jLo, jHi)
        END IF
        CALL GADGOP (BD,NAS,'+')
        bStat = GA_Destroy (lg_B)
        DO I=1,NAS
          SDiag=SD(I)+1.0E-15_wp
          BD(I)=BD(I)/SDiag
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,BD,NAS,IDB)
        CALL mma_deallocate(BD)
! Write the transformation matrices after the diagonal values of B.
! FIXME: This should be removed when the transformation matrices are
! stored as disk resident arrays only.
        IF (KING()) THEN
          CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
          call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,TRANS,1)
            WRITE(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
          CALL mma_deallocate(TRANS)
        END IF
        bStat = GA_Destroy (lg_T)
        RETURN
      END IF
      CALL mma_deallocate(SD)

! TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_T,NAS)
        WRITE(u6,'(1X,A,ES21.14)')                                      &
     &   'EIGENVECTOR NORM BEFORE B TRANS: ', FP
      END IF

      CALL PSBMAT_GETMEM ('B',lg_B,NAS)
      CALL PSBMAT_READ ('B',iCase,iSym,lg_B,NAS)

      If ((do_grad.or.nStpGrd.eq.2).and.do_lindep) Then
        !! The original B matrix is needed in the LinDepLag subroutine
        IDB2 = idBoriMat(ISYM,ICASE)
        call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
        call GA_Access (lg_B, iLo, iHi, jLo, jHi, mB, LDB)
        CALL DDAFILE(LUSTD,1,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB2)
        call GA_Release (lg_B, iLo, iHi, jLo, jHi)
      End If

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(u6,'(1X,A,ES21.14)') 'BMAT NORM: ', FP
      END IF

! FIXME: Perform transformation of B using horizontal stripes of B or
! vertical stripes of T to reduce memory usage if necessary as indicated
! by the available memory, which is now scaling as approx. 3*(NAS**2).
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, One,                      &
     &               lg_B, lg_T, Zero, lg_X )
      bStat = GA_Destroy (lg_B)

      CALL GA_CREATE_STRIPED ('H',NIN,NIN,'BMAT',lg_B)
      call GA_DGEMM ('T', 'N', NIN, NIN, NAS, One,                      &
     &               lg_T, lg_X, Zero, lg_B )
      bStat = GA_Destroy (lg_X)
      bStat = GA_Destroy (lg_T)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_B,NIN)
        WRITE(u6,'(1X,A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

      CALL TIMING(CPU1,CPUE,TIO,TIOE)

! Diagonalize the transformed B matrix.
      CALL mma_allocate(EIG,NIN,Label='EIG')
      EIG(:)=Zero
      IF(BSPECT.NE.'YES')  THEN
! Use diagonal approxim., if allowed.
!        call GA_Fill (lg_V, Zero)
        call GA_Zero (lg_V)
! FIXME: this original code seemed wrong, using uninitialized SD?
!       IDIAG=1
!       DO I=1,NIN
!         EIG(I)=B(IDIAG)/SD
!         IDIAG=IDIAG+1+NIN-I
!       END DO
        WRITE(u6,*) 'GLOB_SBDIAG: option not implemented'
        call AbEnd()
      ELSE
#ifdef _SCALAPACK_
        CALL GA_CREATE_STRIPED ('H',NIN,NIN,'VMAT',lg_V)
        CALL GA_PDSYEVX_ (lg_B, lg_V, EIG, 0)
        bStat = GA_Destroy (lg_B)
#else
        IF (myRank.EQ.0) THEN
          CALL mma_allocate(VEC,NIN**2,Label='VEC')
          CALL GA_Get (lg_B, 1, NIN, 1, NIN, VEC, NIN)
        END IF
        bSTAT = GA_Destroy (lg_B)
        IF (myRank.EQ.0) THEN
          call dsyev_('V','L',NIN,VEC,NIN,EIG,WGRONK,-1,INFO)
          NSCRATCH=INT(WGRONK(1))
          CALL mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
          call dsyev_('V','L',NIN,VEC,NIN,EIG,                          &
     &               SCRATCH,NSCRATCH,INFO)
          CALL mma_deallocate(SCRATCH)
        END IF
        call GA_Sync()
        CALL GA_CREATE_STRIPED ('H',NIN,NIN,'VMAT',lg_V)
        IF (myRank.EQ.0) THEN
          CALL GA_Put (lg_V, 1, NIN, 1, NIN, VEC, NIN)
          CALL mma_deallocate(VEC)
        END IF
        CALL GADGOP(EIG,NIN,'+')
#endif
      END IF

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NIN,EIG,1)
        WRITE(u6,'(1X,A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

! The eigenvalues are written back at same position as the
! original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,EIG,NIN,IDB)
      CALL mma_deallocate(EIG)

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1

! Finally, we must form the composite transformation matrix: T(NAS,NIN)
! matrix on disk * V(NIN*NIN) matrix in core.  FIXME: for now, asume
! there is enough memory for the full transformation, scaling as
! approx. 3*(NAS**2).  Should be determined by the available memory.
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      CALL PSBMAT_READ ('T',iCase,iSym,lg_X,NAS*NIN)
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
      call GA_DGEMM ('N', 'N', NAS, NIN, NIN, One,                      &
     &               lg_X, lg_V, Zero, lg_T )
      bStat = GA_Destroy (lg_X)
      bStat = GA_Destroy (lg_V)

! Write the composite transformation matrix to disk.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

! Additonally, compute S*T and store in on disk for later use by the RHS
! vector utitlities
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)

      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'STMAT',lg_ST)
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, One,                      &
     &               lg_S, lg_T, Zero, lg_ST )
      bStat = GA_Destroy (lg_S)

      CALL PSBMAT_WRITE ('M',iCase,iSym,lg_ST,NAS*NIN)
      bStat = GA_Destroy (lg_ST)

! For now, also keep the transformation matrix on disk as a
! replicate array.  FIXME: Should be removed later.
      IF (KING()) THEN
        CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
        call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
        dTRANS=dNRM2_(NAS*NIN,TRANS,1)
        IF (iPrGlb.GE.INSANE) THEN
          WRITE(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
        END IF
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
        CALL mma_deallocate(TRANS)
      END IF

      call ga_sync()
      bStat = GA_Destroy (lg_T)

#include "macros.fh"
      unused_var(bStat)

      END SUBROUTINE SBDIAG_MPP

#elif defined (NAGFOR)
      ! Some compilers do not like empty files
      subroutine empty_sbmiag_mpp()
      end subroutine empty_sbmiag_mpp
#endif
