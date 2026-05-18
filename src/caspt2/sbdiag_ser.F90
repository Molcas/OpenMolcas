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
      SUBROUTINE SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
      use constants, only: Zero, One
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_grad, do_lindep, nStpGrd, LUSTD,      &
     &                         idBoriMat
      use caspt2_global, only: LUSOLV, LUSBT
      use PrintLevel, only: INSANE
      use EQSOLV, only: IDTMAT, IDBMAT, IDSMAT, IDSTMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: BMatrix, BSpect, BTrans, IfDOrtho,       &
     &                         ThrShn, ThrShs, nASup, nISup, Cases,     &
     &                         nInDep, nG3
      use definitions, only: wp, iwp, ItoB, u6
      IMPLICIT None

      integer(kind=iwp), Intent(in):: iSym, iCase
      real(kind=wp), Intent(out):: CondNr, CPU

! For fooling some compilers:
      REAL(kind=wp) WGRONK(2)

      REAL(kind=wp), allocatable:: S(:), SD(:), SCA(:)
      REAL(kind=wp), allocatable:: VEC(:), EIG(:), SCRATCH(:)
      REAL(kind=wp), allocatable:: B(:), BD(:), BX(:), XBX(:)
      REAL(kind=wp), allocatable:: TRANS(:), AUX(:), ST(:)

      REAL(kind=wp) :: CPU1, CPU2, CPUE, EVAL, FACT, FP, SDIAG, SZ,     &
     &                 SZMAX, SZMIN, TIO, TIOE
      REAL(kind=wp), external :: DNRM2_
      integer(kind=iwp) :: I, IDB, IDB2, IDIAG, IDS, IDST, IDT, IDTMP,  &
     &                     IDTMP0, IJ, INFO, iPad, J, KEND, KSTA,       &
     &                     LTRANS1, LVNEW, LVSTA, NAS, NAUX, NB, NBNEW, &
     &                     NCOEF, NCOL, NIN, NIS, NS, NSCRATCH, JOFF

! On entry, the file LUSBT contains overlap matrices SMAT at disk
! addresses IDSMAT(ISYM,ICASE), ISYM=1..NSYM, ICASE=1..11, and
! similarly BMAT matrices at IDBMAT(ISYM,ICASE).
! These matrices are stored in a triangular format.
! The rectangular matrices TRANS are computed, such that
!    Sum(J) of BMAT(I,J)*TRANS(J,MU) =
!       Sum(J) of  SMAT(I,J)*TRANS(J,MU)*BDIAG(MU)
! where I,J are in 1..NASUP(ISYM,ICASE)
!        MU is in 1..NINDEP(ISYM,ICASE)
! NINDEP is the numerically effective rank of SMAT, and the
! columns of TRANS are orthonormal:
!    Sum(I,J) of  SMAT(I,J)*TRANS(J,MU)*TRANS(J,NU) = Kron(MU,NU)
! BMAT is destroyed, and is overwritten by BDIAG(MU) and
! TRANS(I,MU), with addresses IDBMAT(ISYM,ICASE) and
! IDTMAT(ISYM,ICASE). Enough file space was thus originally
! reserved on LUSBT to allow this overlay.
! LUSOLV is assumed not to be in use yet, so we use it
! for temporary storage.

      SDiag = Zero ! dummy initialize

      CPU=Zero
      CONDNR=Zero
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS

      IF(NCOEF.EQ.0) RETURN

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(u6,'("DEBUG> ",A12,A7,I2,A2,A6,A2,A5,I1)')                &
     &  'SBDIAG_SER: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM
      END IF

      IDTMP0 = 0
      If (do_grad.or.nStpGrd.EQ.2) Then
        !! correct?
        iPad = ItoB - Mod(6*NG3,ItoB)
        IDTMP0=6*NG3+iPad
      End If

      NS=(NAS*(NAS+1))/2
      CALL mma_allocate(S,NS,Label='S')
      IDS=IDSMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,2,S,NS,IDS)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NS,S,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ', FP
      END IF

! For some purposes, we need to save the diagonal elements:
      IF(BMATRIX.EQ.'YES') THEN
        IF(BTRANS.NE.'YES') THEN
          CALL mma_allocate(SD,NAS,Label='SD')
          IDIAG=0
          DO I=1,NAS
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
        END IF
      END IF

! FIRST, FIND NIN ORTHONORMAL VECTORS BY SCALED SYMMETRIC ON.
! Addition, for the scaled symmetric ON: the S matrix is scaled
! to make the diagonal elements close to 1.
! Extremely small values give scale factor exactly zero.
      CALL mma_allocate(SCA,NAS,Label='SCA')
      IDIAG=0
      DO I=1,NAS
        IDIAG=IDIAG+I
        SDiag=S(IDIAG)
        If (IFDORTHO) then
          SCA(I)=One
        Else
          IF(SDiag.GT.THRSHN) THEN
! Small variations of the scale factor were beneficial
            SCA(I)=(One+DBLE(I)*3.0E-6_wp)/SQRT(SDiag)
          ELSE
            SCA(I)=Zero
          END IF
        End If
      END DO
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          S(IJ)=S(IJ)*SCA(I)*SCA(J)
        END DO
      END DO
! End of addition.
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NS,S,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ',  &
     &        FP
      END IF

! DIAGONALIZE THE SCALED S MATRIX:
      CALL mma_allocate(VEC,NAS**2,Label='VEC')
      CALL mma_allocate(EIG,NAS,Label='EIG')

      CALL TIMING(CPU1,CPUE,TIO,TIOE)
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          VEC(NAS*(J-1)+I)=S(IJ)
        END DO
      END DO
      INFO=0
      call dsyev_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
      NSCRATCH=INT(WGRONK(1))
      CALL mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
      call dsyev_('V','U',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
      CALL mma_deallocate(SCRATCH)
      CALL mma_deallocate(S)

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1
      ! fingerprint eigenvalues
      if (iprglb >= insane) then
        fp = dnrm2_(nas,eig,1)
        write(u6,'("DEBUG> ",A,ES21.14)') 'Smat eigval norm: ', fp
      end if

! Form orthonormal vectors by scaling eigenvectors
      NIN=0
      DO I=1,NAS
        EVAL=EIG(I)
        IF(EVAL.LT.THRSHS) CYCLE
        FACT=One/SQRT(EVAL)
        NIN=NIN+1
        LVSTA=1+NAS*(I-1)
        IF(NIN.EQ.I) THEN
          CALL DSCAL_(NAS,FACT,VEC(LVSTA:),1)
        ELSE
          LVNEW=1+NAS*(NIN-1)
          CALL DYAX(NAS,FACT,VEC(LVSTA:),1,VEC(LVNEW:),1)
        END IF
      END DO
      NINDEP(ISYM,ICASE)=NIN
      CALL mma_deallocate(EIG)
! Addition, for the scaled symmetric ON.
      DO I=1,NAS
        CALL DSCAL_(NIN,SCA(I),VEC(I:),NAS)
      END DO

      CALL mma_deallocate(SCA)
! The condition number, after scaling, disregarding linear dep.
      IF(NIN.GE.2) THEN
        SZMIN=1.0E99_wp
        SZMAX=Zero
        DO I=1,NIN
          SZ=DNRM2_(NAS,VEC(1+NAS*(I-1):),1)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=(SZMAX/SZMIN)**2
      END IF
! End of addition.
      IF(NIN.EQ.0) THEN
        CALL mma_deallocate(VEC)
        RETURN
      END IF

      IF(BMATRIX.EQ.'NO') THEN
! In some calculations, we do not use B matrices.
! Just write the transformation matrix and branch out:
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
        CALL mma_deallocate(VEC)
        IF (IPRGLB.GE.INSANE) THEN
          WRITE(u6,'("DEBUG> ",A)') 'SBDIAG: skip B matrix'
        END IF
        RETURN
      ELSE IF (BTRANS.NE.'YES') THEN
! In other calculations, B matrix is used, but not transformed.
! We may need the diagonal active energies: the diagonal values of
! B divided by the diagonal values of S. These are placed where
! the eigenvalues would go in ordinary CASPT2.
! NOTE: On LUSBT, the transformation matrices partly overwrite
! and destroy the B matrices. The diagonal elements of B must be
! extracted before the transformation matrix is written.
        CALL mma_allocate(BD,NAS,Label='BD')
        NB=(NAS*(NAS+1))/2
        CALL mma_allocate(B,NB,Label='B')
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,B,NB,IDB)
        IDIAG=0
        DO I=1,NAS
          IDIAG=IDIAG+I
          BD(I)=B(IDIAG)
        END DO
        CALL mma_deallocate(B)
! Now, the transformation matrix can be written out.
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
        CALL mma_deallocate(VEC)
        DO I=1,NAS
          SDiag=SD(I)+1.0e-15_wp
          BD(I)=BD(I)/SDiag
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,BD,NAS,IDB)
        CALL mma_deallocate(SD)
        CALL mma_deallocate(BD)
        IF (IPRGLB.GE.INSANE) THEN
          WRITE(u6,'("DEBUG> ",A)')                                     &
     &         'SBDIAG: skip B matrix transformation'
          WRITE(u6,'("DEBUG> ",A)')'        but keep B_ii/S_ii values'
        END IF
        RETURN
      END IF

! TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
! USE LUSOLV AS TEMPORARY STORAGE. WE MAY NEED SECTIONING.
! NOTE: SECTIONING MUST BE  PRECISELY THE SAME AS WHEN LATER
! READ BACK (SEE BELOW).
      NAUX=MIN(19,NIN)
      IDTMP=IDTMP0
      CALL DDAFILE(LUSOLV,1,VEC,NAS*NAUX,IDTMP)
      DO KSTA=NAUX+1,NIN,NAUX
        KEND=MIN(KSTA-1+NAUX,NIN)
        NCOL=1+KEND-KSTA
        LVSTA=1+NAS*(KSTA-1)
        CALL DDAFILE(LUSOLV,1,VEC(LVSTA),NAS*NCOL,IDTMP)
      END DO
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS**2,VEC,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)')                               &
     &   'EIGENVECTOR NORM BEFORE B TRANS: ', FP
      END IF

! TRANSFORM B. NEW B WILL OVERWRITE AND DESTROY VEC
      IDB=IDBMAT(ISYM,ICASE)
      NB=NS
      CALL mma_allocate(B,NB,Label='B')
      CALL DDAFILE(LUSBT,2,B,NB,IDB)
      If ((do_grad.or.nStpGrd.eq.2).and.do_lindep) Then
        !! The original B matrix is needed in the LinDepLag subroutine
        IDB2 = idBoriMat(ISYM,ICASE)
        CALL DDAFILE(LUSTD,1,B,NB,IDB2)
      End If
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NB,B,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM: ', FP
      END IF

      CALL mma_allocate(BX,NAS,Label='BX')
      CALL mma_allocate(XBX,NAS,Label='XBX')
      DO J=NIN,1,-1
        LVSTA=1+NAS*(J-1)
        BX(:)=Zero
#ifdef _CRAY_C90_
        CALL SSPMV('U',NAS,One,B,VEC(LVSTA),1,                          &
     &                           One,BX,1)
#else
!        CALL DSLMX(NAS,One,B,VEC(LVSTA),1,
!     &                                   BX,1)
        CALL DSPMV_('U',NAS,One,B,VEC(LVSTA),1,                         &
     &                           One,BX,1)
#endif
! BX: B * Vector number J.
        CALL DCOPY_(J,[Zero],0,XBX,1)
        CALL DGEMM_('T','N',                                            &
     &              J,1,NAS,                                            &
     &              One,VEC,NAS,                                        &
     &              BX,NAS,                                             &
     &              Zero,XBX,J)
! XBX CONTAINS NOW THE UPPERTRIANGULAR
! ELEMENTS OF THE J-th COLUMN OF TRANSFORMED B MATRIX.
        CALL DCOPY_(J,XBX,1,VEC(LVSTA),1)
      END DO
      CALL mma_deallocate(BX)
      CALL mma_deallocate(XBX)
      CALL mma_deallocate(B)
! VEC HAS NOW BEEN DESTROYED (OVERWRITTEN BY NEW B).
! COPY TO TRIANGULAR STORAGE.
      NBNEW=(NIN*(NIN+1))/2
      CALL mma_allocate(B,NBNEW,Label='B')
      DO J=1,NIN
        JOFF=(J*(J-1))/2
        CALL DCOPY_(J,VEC(1+NAS*(J-1):),1,B(1+JOFF:),1)
      END DO
      CALL mma_deallocate(VEC)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NBNEW,B,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

! DIAGONALIZE THE TRANSFORMED B MATRIX.
      CALL mma_allocate(EIG,NIN,Label='EIG')
      CALL mma_allocate(VEC,NIN**2,Label='VEC')
      CALL TIMING(CPU1,CPUE,TIO,TIOE)
! - Alt 0: Use diagonal approxim., if allowed:
      IF(BSPECT.NE.'YES')  THEN
        IDIAG=0
        Write (6,*) 'Sbdiag: THis code does not make sense!'
        Write (6,*) '        SDiag is not properly defined!'
        Call Abend()
        DO I=1,NIN
          IDIAG=IDIAG+I
          EIG(I)=B(IDIAG)/SDiag
        END DO
      ELSE
        IJ=0
        DO J=1,NIN
          DO I=1,J
            IJ=IJ+1
            VEC(NIN*(J-1)+I)=B(IJ)
          END DO
        END DO
        CALL DSYEV_('V','U',NIN,VEC,NIN,EIG,WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        CALL mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
        CALL DSYEV_('V','U',NIN,VEC,NIN,EIG,SCRATCH,NSCRATCH,INFO)
        CALL mma_deallocate(SCRATCH)
        CALL mma_deallocate(B)
      END IF
      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NIN,EIG,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

! The eigenvalues are written back at same position as the
! original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,EIG,NIN,IDB)
      CALL mma_deallocate(EIG)

! Finally, we must form the composite transformation matrix,
!  = (NAS*NIN matrix on disk) * (NIN*NIN matrix in core).
! Assume enough space since we got rid of S/B matrices.
! Specifically, assume we have enough space for the two
! full matrices, plus an additional 19 columns of results.
      NAUX=MIN(19,NIN)
      CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
      CALL mma_allocate(AUX  ,NAS*NAUX,Label='AUX')
      IDTMP=IDTMP0
      CALL DDAFILE(LUSOLV,2,AUX,NAS*NAUX,IDTMP)
      IF(BTRANS.EQ.'YES') THEN
        CALL DGEMM_('N','N',                                            &
     &              NAS,NIN,NAUX,                                       &
     &              One,AUX,NAS,                                        &
     &              VEC,NIN,                                            &
     &              Zero,TRANS,NAS)
      ELSE
        CALL DCOPY_(NAS*NAUX,AUX,1,TRANS,1)
      END IF
      DO KSTA=NAUX+1,NIN,NAUX
        KEND=MIN(KSTA-1+NAUX,NIN)
        NCOL=1+KEND-KSTA
        CALL DDAFILE(LUSOLV,2,AUX,NAS*NCOL,IDTMP)
        IF(BTRANS.EQ.'YES') THEN
          CALL DGEMM_('N','N',NAS,NIN,NCOL,One,                         &
     &              AUX,NAS,VEC(KSTA),NIN,                              &
     &              One,TRANS,NAS)
        ELSE
          LTRANS1=1+NAS*(KSTA-1)
          CALL DCOPY_(NAS*NCOL,AUX,1,TRANS(LTRANS1),1)
        END IF
      END DO
      CALL mma_deallocate(AUX)
      CALL mma_deallocate(VEC)
      IDT=IDTMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS*NIN,TRANS,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', FP
      END IF

!-SVC: compute S*T and store on disk for later use by RHS vector
!      utilities.
      NS=(NAS*(NAS+1))/2
      CALL mma_allocate(S,NS,Label='S')
      IDS=IDSMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,2,S,NS,IDS)
      CALL mma_allocate(ST,NAS*NIN,Label='ST')
      ST(:)=Zero
      CALL TRIMUL(NAS,NIN,One,S,TRANS,NAS,ST,NAS)
      CALL mma_deallocate(S)
      CALL mma_deallocate(TRANS)
      IDST=IDSTMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,ST,NAS*NIN,IDST)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS*NIN,ST,1)
        WRITE(u6,'("DEBUG> ",A,ES21.14)') 'STMAT NORM: ', FP
      END IF
      CALL mma_deallocate(ST)

      END SUBROUTINE SBDIAG_SER
