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
      SUBROUTINE SBDIAG()
      use definitions, only: iwp, wp
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: USUAL, VERBOSE
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_module, only: nSym, ThrShn, ThrShs, Cases, nASup,
     &  nISup, nInDep
      IMPLICIT None

      real(kind=wp) CondNr, CPU
      integer(kind=iwp) iCase, iPar0, iPar1, iSym


      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' Find transformation matrices to eigenbasis'//
     &     ' of block-diagonal part of H0.'
        WRITE(6,*)' Eliminate linear dependency. Thresholds for:'
        WRITE(6,'(A,G12.4)')'   Initial squared norm  :',THRSHN
        WRITE(6,'(A,G12.4)')'   Eigenvalue of scaled S:',THRSHS
      END IF

C SVC.20100904: there are now two SBDIAG versions: a replicate
C subroutine expecting replicate S and B matrices (in upper triangular
C column-wise storage) and performing transformations on each process,
C and a global array subroutine expecting S and B matrices in global
C arrays (full column-wise storage) spread over processes.  The latter
C is currently used only for cases 1 (A) and 4 (C), for which global
C array mksmat and mkbmat routines have been implemented, as these can
C grow very big with increasing size of the active space.  For now, we
C still use the replicate routines for the other cases as they have more
C modest array sizes.

      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' Condition numbers are computed after diagonal'//
     &     ' scaling and after removal of'
        WRITE(6,*)' linear dependency. Resulting sizes, condition'//
     &     ' numbers, and times:'
        WRITE(6,'(3X,A10,4A12,A9)')
     &     'CASE(SYM)','NASUP','NISUP','NINDEP','COND NR','CPU (s)'
      ENDIF

      DO iCASE=1,11
        DO ISYM=1,NSYM
#ifdef _MOLCAS_MPP_
            IF (IS_REAL_PAR() .AND.
     &          (ICASE.EQ.1.OR.ICASE.EQ.4)) THEN
              CALL SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
            ELSE
#endif
              CALL SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
#ifdef _MOLCAS_MPP_
            END IF
#endif
            IF (IPRGLB.GE.VERBOSE) THEN
              WRITE(6,'(3X,A6,A1,I1,A1,1X,3I12,G11.2,I9)')
     &         CASES(ICASE),'(',ISYM,')',
     &         NASUP(ISYM,ICASE),NISUP(ISYM,ICASE),
     &         NINDEP(ISYM,ICASE),CONDNR,NINT(CPU)
            END IF

        END DO
      END DO

C usually print info on the total number of parameters
      IPAR0=0
      IPAR1=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          IPAR0=IPAR0+NASUP(ISYM,ICASE)*NISUP(ISYM,ICASE)
          IPAR1=IPAR1+NINDEP(ISYM,ICASE)*NISUP(ISYM,ICASE)
        END DO
      END DO
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)' Total nr of CASPT2 parameters:'
        WRITE(6,'(a,i12)')'   Before reduction:',IPAR0
        WRITE(6,'(a,i12)')'   After  reduction:',IPAR1
      ENDIF

      END SUBROUTINE SBDIAG

      SUBROUTINE SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
      use definitions, only: wp, iwp, ItoB
      use constants, only: Zero, One
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_grad, do_lindep, nStpGrd, LUSTD,
     *                         idBoriMat
      use caspt2_global, only: LUSOLV, LUSBT
      use PrintLevel, only: INSANE
      use EQSOLV, only: IDTMAT, IDBMAT, IDSMAT, IDSTMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: BMatrix, BSpect, BTrans, IfDOrtho,
     &                         ThrShn, ThrShs, nASup, nISup, Cases,
     &                         nInDep

      use pt2_guga, only: nG3
      IMPLICIT None

      integer(kind=iwp), Intent(in):: iSym, iCase
      real(kind=wp), Intent(out):: CondNr, CPU

* For fooling some compilers:
      REAL(kind=wp) WGRONK(2)

      REAL(kind=wp), allocatable:: S(:), SD(:), SCA(:)
      REAL(kind=wp), allocatable:: VEC(:), EIG(:), SCRATCH(:)
      REAL(kind=wp), allocatable:: B(:), BD(:), BX(:), XBX(:)
      REAL(kind=wp), allocatable:: TRANS(:), AUX(:), ST(:)

      REAL(kind=wp) :: CPU1, CPU2, CPUE, EVAL, FACT, FP, SDIAG, SZ,
     &                 SZMAX, SZMIN, TIO, TIOE
      REAL(kind=wp), external :: DNRM2_
      integer(kind=iwp) :: I, IDB, IDB2, IDIAG, IDS, IDST, IDT, IDTMP,
     &                     IDTMP0, IJ, INFO, iPad, J, KEND, KSTA,
     &                     LTRANS1, LVNEW, LVSTA, NAS, NAUX, NB, NBNEW,
     &                     NCOEF, NCOL, NIN, NIS, NS, NSCRATCH, JOFF

C On entry, the file LUSBT contains overlap matrices SMAT at disk
C addresses IDSMAT(ISYM,ICASE), ISYM=1..NSYM, ICASE=1..11, and
C similarly BMAT matrices at IDBMAT(ISYM,ICASE).
C These matrices are stored in a triangular format.
C The rectangular matrices TRANS are computed, such that
C    Sum(J) of BMAT(I,J)*TRANS(J,MU) =
C       Sum(J) of  SMAT(I,J)*TRANS(J,MU)*BDIAG(MU)
C where I,J are in 1..NASUP(ISYM,ICASE)
C        MU is in 1..NINDEP(ISYM,ICASE)
C NINDEP is the numerically effective rank of SMAT, and the
C columns of TRANS are orthonormal:
C    Sum(I,J) of  SMAT(I,J)*TRANS(J,MU)*TRANS(J,NU) = Kron(MU,NU)
C BMAT is destroyed, and is overwritten by BDIAG(MU) and
C TRANS(I,MU), with addresses IDBMAT(ISYM,ICASE) and
C IDTMAT(ISYM,ICASE). Enough file space was thus originally
C reserved on LUSBT to allow this overlay.
C LUSOLV is assumed not to be in use yet, so we use it
C for temporary storage.

      SDiag = Zero ! dummy initialize

      CPU=Zero
      CONDNR=Zero
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS

      IF(NCOEF.EQ.0) RETURN

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(6,'("DEBUG> ",A12,A7,I2,A2,A6,A2,A5,I1)')
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ', FP
      END IF

C For some purposes, we need to save the diagonal elements:
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

C FIRST, FIND NIN ORTHONORMAL VECTORS BY SCALED SYMMETRIC ON.
C Addition, for the scaled symmetric ON: the S matrix is scaled
C to make the diagonal elements close to 1.
C Extremely small values give scale factor exactly zero.
      CALL mma_allocate(SCA,NAS,Label='SCA')
      IDIAG=0
      DO I=1,NAS
        IDIAG=IDIAG+I
        SDiag=S(IDIAG)
        If (IFDORTHO) then
          SCA(I)=One
        Else
          IF(SDiag.GT.THRSHN) THEN
* Small variations of the scale factor were beneficial
            SCA(I)=(One+DBLE(I)*3.0D-6)/SQRT(SDiag)
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
C End of addition.
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NS,S,1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ', FP
      END IF

C DIAGONALIZE THE SCALED S MATRIX:
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
        write(6,'("DEBUG> ",A,ES21.14)') 'Smat eigval norm: ', fp
      end if

C Form orthonormal vectors by scaling eigenvectors
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
C Addition, for the scaled symmetric ON.
      DO I=1,NAS
        CALL DSCAL_(NIN,SCA(I),VEC(I:),NAS)
      END DO

      CALL mma_deallocate(SCA)
C The condition number, after scaling, disregarding linear dep.
      IF(NIN.GE.2) THEN
        SZMIN=1.0D99
        SZMAX=Zero
        DO I=1,NIN
          SZ=DNRM2_(NAS,VEC(1+NAS*(I-1):),1)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=(SZMAX/SZMIN)**2
      END IF
C End of addition.
      IF(NIN.EQ.0) THEN
        CALL mma_deallocate(VEC)
        RETURN
      END IF

      IF(BMATRIX.EQ.'NO') THEN
C In some calculations, we do not use B matrices.
C Just write the transformation matrix and branch out:
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
        CALL mma_deallocate(VEC)
        IF (IPRGLB.GE.INSANE) THEN
          WRITE(6,'("DEBUG> ",A)') 'SBDIAG: skip B matrix'
        END IF
        RETURN
      ELSE IF (BTRANS.NE.'YES') THEN
C In other calculations, B matrix is used, but not transformed.
C We may need the diagonal active energies: the diagonal values of
C B divided by the diagonal values of S. These are placed where
C the eigenvalues would go in ordinary CASPT2.
C NOTE: On LUSBT, the transformation matrices partly overwrite
C and destroy the B matrices. The diagonal elements of B must be
C extracted before the transformation matrix is written.
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
C Now, the transformation matrix can be written out.
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
        CALL mma_deallocate(VEC)
        DO I=1,NAS
          SDiag=SD(I)+1.0d-15
          BD(I)=BD(I)/SDiag
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,BD,NAS,IDB)
        CALL mma_deallocate(SD)
        CALL mma_deallocate(BD)
        IF (IPRGLB.GE.INSANE) THEN
          WRITE(6,'("DEBUG> ",A)')'SBDIAG: skip B matrix transformation'
          WRITE(6,'("DEBUG> ",A)')'        but keep B_ii/S_ii values'
        END IF
        RETURN
      END IF

C TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
C USE LUSOLV AS TEMPORARY STORAGE. WE MAY NEED SECTIONING.
C NOTE: SECTIONING MUST BE  PRECISELY THE SAME AS WHEN LATER
C READ BACK (SEE BELOW).
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
        WRITE(6,'("DEBUG> ",A,ES21.14)')
     &   'EIGENVECTOR NORM BEFORE B TRANS: ', FP
      END IF

C TRANSFORM B. NEW B WILL OVERWRITE AND DESTROY VEC
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM: ', FP
      END IF

      CALL mma_allocate(BX,NAS,Label='BX')
      CALL mma_allocate(XBX,NAS,Label='XBX')
      DO J=NIN,1,-1
        LVSTA=1+NAS*(J-1)
        BX(:)=Zero
#ifdef _CRAY_C90_
        CALL SSPMV('U',NAS,One,B,VEC(LVSTA),1,
     &                           One,BX,1)
#else
*        CALL DSLMX(NAS,One,B,VEC(LVSTA),1,
*     &                                   BX,1)
        CALL DSPMV_('U',NAS,One,B,VEC(LVSTA),1,
     &                           One,BX,1)
#endif
C BX: B * Vector number J.
        CALL DCOPY_(J,[Zero],0,XBX,1)
        CALL DGEMM_('T','N',
     &              J,1,NAS,
     &              One,VEC,NAS,
     &              BX,NAS,
     &              Zero,XBX,J)
C XBX CONTAINS NOW THE UPPERTRIANGULAR
C ELEMENTS OF THE J-th COLUMN OF TRANSFORMED B MATRIX.
        CALL DCOPY_(J,XBX,1,VEC(LVSTA),1)
      END DO
      CALL mma_deallocate(BX)
      CALL mma_deallocate(XBX)
      CALL mma_deallocate(B)
C VEC HAS NOW BEEN DESTROYED (OVERWRITTEN BY NEW B).
C COPY TO TRIANGULAR STORAGE.
      NBNEW=(NIN*(NIN+1))/2
      CALL mma_allocate(B,NBNEW,Label='B')
      DO J=1,NIN
        JOFF=(J*(J-1))/2
        CALL DCOPY_(J,VEC(1+NAS*(J-1):),1,B(1+JOFF:),1)
      END DO
      CALL mma_deallocate(VEC)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NBNEW,B,1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

C DIAGONALIZE THE TRANSFORMED B MATRIX.
      CALL mma_allocate(EIG,NIN,Label='EIG')
      CALL mma_allocate(VEC,NIN**2,Label='VEC')
      CALL TIMING(CPU1,CPUE,TIO,TIOE)
C - Alt 0: Use diagonal approxim., if allowed:
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

C The eigenvalues are written back at same position as the
C original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,EIG,NIN,IDB)
      CALL mma_deallocate(EIG)

C Finally, we must form the composite transformation matrix,
C  = (NAS*NIN matrix on disk) * (NIN*NIN matrix in core).
C Assume enough space since we got rid of S/B matrices.
C Specifically, assume we have enough space for the two
C full matrices, plus an additional 19 columns of results.
      NAUX=MIN(19,NIN)
      CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
      CALL mma_allocate(AUX  ,NAS*NAUX,Label='AUX')
      IDTMP=IDTMP0
      CALL DDAFILE(LUSOLV,2,AUX,NAS*NAUX,IDTMP)
      IF(BTRANS.EQ.'YES') THEN
        CALL DGEMM_('N','N',
     &              NAS,NIN,NAUX,
     &              One,AUX,NAS,
     &              VEC,NIN,
     &              Zero,TRANS,NAS)
      ELSE
        CALL DCOPY_(NAS*NAUX,AUX,1,TRANS,1)
      END IF
      DO KSTA=NAUX+1,NIN,NAUX
        KEND=MIN(KSTA-1+NAUX,NIN)
        NCOL=1+KEND-KSTA
        CALL DDAFILE(LUSOLV,2,AUX,NAS*NCOL,IDTMP)
        IF(BTRANS.EQ.'YES') THEN
          CALL DGEMM_('N','N',NAS,NIN,NCOL,One,
     &              AUX,NAS,VEC(KSTA),NIN,
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', FP
      END IF

C-SVC: compute S*T and store on disk for later use by RHS vector
C      utilities.
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'STMAT NORM: ', FP
      END IF
      CALL mma_deallocate(ST)

      END SUBROUTINE SBDIAG_SER

C SVC2010: The following subroutine computes the transformation
C matrices using global arrays rather than replicate arrays.  Currently,
C the maximum memory scales approximately as 3*(NAS**2), which is twice
C what is needed in the replicate routines.  It can be reduced to
C 2*(NAS**2), if the matrix multiplications are rewritten to occur in
C batch mode.  However, unlike in the replicate routine, this amount is
C divided over processors.
#ifdef _MOLCAS_MPP_
      SUBROUTINE SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
      use definitions, only: iwp, wp
      use Constants, only: Zero, One
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: INSANE
      USE Para_Info, ONLY: King
      use caspt2_global, only: LUSBT
      use caspt2_global, only: do_grad, do_lindep, nStpGrd, LUSTD,
     &                         idBoriMat
      use EQSOLV, only: IDSMAT, IDTMAT, IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nASup, nISup, Cases, IfDOrtho, ThrShn,
     &                         ThrShs, nInDep, BMATRIX, BTRANS, BSPECT
      IMPLICIT None

      integer(kind=iwp), intent(in):: iSym, iCase
      real(kind=wp), intent(out):: CondNr, CPU
C-SVC20100902: global arrays header files
#include "global.fh"
#include "mafdecls.fh"
#ifndef _SCALAPACK_
      real(kind=wp) WGRONK(2)
#endif
      LOGICAL(kind=iwp) bSTAT
      CHARACTER(LEN=2) cSYM,cCASE
      real(kind=wp), allocatable:: COL(:), TMP(:), SD(:), SCA(:),
     &                             EIG(:), VEC(:), SCRATCH(:), COND(:),
     &                             TRANS(:), BD(:)
      integer(kind=iwp) NAS, NIS, NCOEF, lg_S, NCOL, NTMP, IOFF, J, IDS,
     &                  MyRank, iLo, iHi, jLo, jHi, ISTA, IEND, MS, LDS,
     &                  I, lg_V, NIN, mV, LDV, lg_T, IDT, lg_B, mB, lDB,
     &                  IDB, IDB2, lg_X, lg_ST, Info, NSCRATCH
      real(kind=wp) FP, SDiag, SZMIN, SZMAX, SZ, dTrans
      real(kind=wp) CPU1, CPUE, TIO, TIOE, CPU2
      real(kind=wp), External:: PSBMAT_FPRINT, DNRM2_

C On entry, the DRA metafiles contain the matrices S and B for cases A
C (iCASE=1) en C (iCASE=4).  These symmetric matrices are stored on disk
C in full square format.  The rectangular matrices T are computed, such
C that Sum_(J) [B(I,J)*T(J,MU)] = Sum_(J) [S(I,J)*T(J,MU)*BD(MU)], where
C I,J are in the range (1,NASUP(ISYM,ICASE)) and MU is in the range
C (1,NINDEP(ISYM,ICASE)).  NINDEP is the numerically effective rank of
C S, and the columns of T are orthonormal: Sum_(I,J)
C [S(I,J)*T(J,MU)*T(J,NU)] = Kron(MU,NU).  B is destroyed, and is
C overwritten by BD(MU) and T(I,MU), which is stored in a DRA metafile.

C Initialize the DRA I/O subsystem with default values.

      IF (iCASE.NE.1.AND.iCASE.NE.4) THEN
        WRITE(6,*) 'Invalid CASE number used for global SBDIAG, Abort'
        CALL AbEnd()
      END IF

      write(unit=cCase, fmt='(I2.2)') iCase
      write(unit=cSYM, fmt='(I2.2)') iSYM
C Start a long loop over irreps:

      CPU=Zero
      CONDNR=Zero
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS
      IF(NCOEF.EQ.0) RETURN

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(6,'("DEBUG> ",A12,A5,I2,A2,A6,A2,A5,I1)')
     &  'SBDIAG_MPP: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM
      END IF

C Allocate memory for the S matrix and read it from disk:
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)
      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_S,NAS)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ', FP
      END IF

C The S matrices are needed later on by non-global routines.  Take the
C oportunity to save them to LUSBT here.  FIXME: Should be removed once
C full parallelization of use of S matrices is achieved.
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

C Save the diagonal elements from the S matrix for easy access later on.
C FIXME: nicer way to do this?
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
      CALL GADSUM (SD,NAS)

C Calculate the scaling factors and store them in array SCA.
      CALL mma_allocate(SCA,NAS,Label='SCA')
      DO I=1,NAS
        SDiag=SD(I)
        IF (IFDORTHO) THEN
          SCA(I)=One
        ELSE
          IF(SDiag.GT.THRSHN) THEN
            SCA(I)=(One+DBLE(I)*3.0D-6)/SQRT(SDiag)
          ELSE
            SCA(I)=Zero
          END IF
        END IF
      END DO

C Scale the elements S(I,J) with the factor SCA(I)*SCA(J).
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
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ', FP
      END IF

      CALL mma_allocate(EIG,NAS,Label='EIG')
      EIG(:)=Zero
C Diagonalize the global array S.  Some (old) reports about parallel
C performance recommend PDSYEVX or PDSYEVR as fastest methods if
C eigenvectors are needed (FIXME: should time this).  For the linear
C dependence removal, split eigenvectors in horizontal stripes so that
C each processor has a row window of all column vectors
#ifdef _SCALAPACK_
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      CALL GA_PDSYEVX_ (lg_S, lg_V, EIG, 0)
      bSTAT = GA_Destroy (lg_S)
#else
C here for the non-ScaLAPACK version: copy matrix to master process,
C diagonalize using the serial DSYEV routine, and copy the resulting
C eigenvectors back to a global array.  Then distribute the eigenvalues.
      IF (myRank.EQ.0) THEN
        CALL mma_allocate(VEC,NAS**2,Label='VEC')
        CALL GA_Get (lg_S, 1, NAS, 1, NAS, VEC, NAS)
      END IF
      bSTAT = GA_Destroy (lg_S)
      IF (myRank.EQ.0) THEN
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        CALL mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,
     &              SCRATCH,NSCRATCH,INFO)
        CALL mma_deallocate(SCRATCH)
      END IF
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      IF (myRank.EQ.0) THEN
        CALL GA_Put (lg_V, 1, NAS, 1, NAS, VEC, NAS)
        CALL mma_deallocate(VEC)
      END IF
      CALL GADSUM(EIG,NAS)
#endif

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS,EIG,1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT EIGENVALUE NORM: ', FP
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
C Form orthonormal transformation vectors by scaling the eigenvectors.
      call GA_Sync()
      myRank = GA_NodeID()
      call GA_Distribution (lg_V, myRank, iLo, iHi, jLo, jHi)
      IF (iLo.NE.0) THEN
        call GA_Access (lg_V, iLo, iHi, jLo, jHi, mV, LDV)
        IF ((jHi-jLo+1).NE.NAS) THEN
          WRITE(6,*) 'SBDIAG_MPP: error in striping of lg_V, ABORT'
          CALL ABEND()
        END IF
        call V_SCALE (EIG,SCA(iLo),DBL_MB(mV),
     &              iHi-iLo+1,jHi-jLo+1,LDV,NIN,COND)
        call GA_Release_Update (lg_V, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()

      CALL mma_deallocate(EIG)
      CALL mma_deallocate(SCA)

C The condition number, after scaling, disregarding linear dep.
C FIXME: adapt to local subroutine for global array lg_V
      IF(NIN.GE.2) THEN
        CALL GADSUM (COND,NIN)
        SZMIN=1.0D99
        SZMAX=Zero
        DO I=1,NIN
          SZ=COND(I)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=SZMAX/SZMIN
      END IF
      CALL mma_deallocate(COND)

C Copy the NIN non-linear dependent eigenvectors to the transformation
C matrix T(NAS,NIN).
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
      call GA_Copy_Patch ('N', lg_V, 1, NAS, 1, NIN,
     &                         lg_T, 1, NAS, 1, NIN)
      bStat = GA_Destroy (lg_V)

      IF(BMATRIX.EQ.'NO      ') THEN
C In some calculations, we do not use B matrices.
C Write the T matrix to disk and exit.  FIXME: This
C should be removed when the transformation matrices are stored as disk
C resident arrays only.
        IF (KING()) THEN
          CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
          call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,TRANS,1)
            WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
          CALL mma_deallocate(TRANS)
        END IF
        bStat = GA_Destroy (lg_T)
        RETURN
      ELSE IF(BTRANS.NE.'YES') THEN
C In other calculations, the B matrix is used but not transformed.  We
C may need the diagonal active energies, i.e. the diagonal values of B
C divided by the diagonal values of S. These are placed where the
C eigenvalues would go in ordinary CASPT2.
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
        CALL GADSUM (BD,NAS)
        bStat = GA_Destroy (lg_B)
        DO I=1,NAS
          SDiag=SD(I)+1.0d-15
          BD(I)=BD(I)/SDiag
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,BD,NAS,IDB)
        CALL mma_deallocate(BD)
C Write the transformation matrices after the diagonal values of B.
C FIXME: This should be removed when the transformation matrices are
C stored as disk resident arrays only.
        IF (KING()) THEN
          CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
          call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,TRANS,1)
            WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
          CALL mma_deallocate(TRANS)
        END IF
        bStat = GA_Destroy (lg_T)
        RETURN
      END IF
      CALL mma_deallocate(SD)

C TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_T,NAS)
        WRITE(6,'(1X,A,ES21.14)')
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
        WRITE(6,'(1X,A,ES21.14)') 'BMAT NORM: ', FP
      END IF

C FIXME: Perform transformation of B using horizontal stripes of B or
C vertical stripes of T to reduce memory usage if necessary as indicated
C by the available memory, which is now scaling as approx. 3*(NAS**2).
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, One,
     &               lg_B, lg_T, Zero, lg_X )
      bStat = GA_Destroy (lg_B)

      CALL GA_CREATE_STRIPED ('H',NIN,NIN,'BMAT',lg_B)
      call GA_DGEMM ('T', 'N', NIN, NIN, NAS, One,
     &               lg_T, lg_X, Zero, lg_B )
      bStat = GA_Destroy (lg_X)
      bStat = GA_Destroy (lg_T)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_B,NIN)
        WRITE(6,'(1X,A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

      CALL TIMING(CPU1,CPUE,TIO,TIOE)

C Diagonalize the transformed B matrix.
      CALL mma_allocate(EIG,NIN,Label='EIG')
      EIG(:)=Zero
      IF(BSPECT.NE.'YES')  THEN
C Use diagonal approxim., if allowed.
C        call GA_Fill (lg_V, 0.0D0)
        call GA_Zero (lg_V)
C FIXME: this original code seemed wrong, using uninitialized SD?
*       IDIAG=1
*       DO I=1,NIN
*         EIG(I)=B(IDIAG)/SD
*         IDIAG=IDIAG+1+NIN-I
*       END DO
        WRITE(6,*) 'GLOB_SBDIAG: option not implemented'
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
          call dsyev_('V','L',NIN,VEC,NIN,EIG,
     &               SCRATCH,NSCRATCH,INFO)
          CALL mma_deallocate(SCRATCH)
        END IF
        call GA_Sync()
        CALL GA_CREATE_STRIPED ('H',NIN,NIN,'VMAT',lg_V)
        IF (myRank.EQ.0) THEN
          CALL GA_Put (lg_V, 1, NIN, 1, NIN, VEC, NIN)
          CALL mma_deallocate(VEC)
        END IF
        CALL GADSUM(EIG,NIN)
#endif
      END IF

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NIN,EIG,1)
        WRITE(6,'(1X,A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

C The eigenvalues are written back at same position as the
C original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,EIG,NIN,IDB)
      CALL mma_deallocate(EIG)

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1

C Finally, we must form the composite transformation matrix: T(NAS,NIN)
C matrix on disk * V(NIN*NIN) matrix in core.  FIXME: for now, asume
C there is enough memory for the full transformation, scaling as
C approx. 3*(NAS**2).  Should be determined by the available memory.
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      CALL PSBMAT_READ ('T',iCase,iSym,lg_X,NAS*NIN)
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
      call GA_DGEMM ('N', 'N', NAS, NIN, NIN, One,
     &               lg_X, lg_V, Zero, lg_T )
      bStat = GA_Destroy (lg_X)
      bStat = GA_Destroy (lg_V)

C Write the composite transformation matrix to disk.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

C Additonally, compute S*T and store in on disk for later use by the RHS
C vector utitlities
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)

      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'STMAT',lg_ST)
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, One,
     &               lg_S, lg_T, Zero, lg_ST )
      bStat = GA_Destroy (lg_S)

      CALL PSBMAT_WRITE ('M',iCase,iSym,lg_ST,NAS*NIN)
      bStat = GA_Destroy (lg_ST)

C For now, also keep the transformation matrix on disk as a
C replicate array.  FIXME: Should be removed later.
      IF (KING()) THEN
        CALL mma_allocate(TRANS,NAS*NIN,Label='TRANS')
        call GA_Get (lg_T, 1, NAS, 1, NIN, TRANS, NAS)
        dTRANS=dNRM2_(NAS*NIN,TRANS,1)
        IF (iPrGlb.GE.INSANE) THEN
          WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
        END IF
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
        CALL mma_deallocate(TRANS)
      END IF

      call ga_sync()
      bStat = GA_Destroy (lg_T)

      END SUBROUTINE SBDIAG_MPP

      SUBROUTINE S_SCALE (NAS,SCA,S,iLo,iHi,jLo,jHi,LDS)
      use definitions, only: iwp, wp

      IMPLICIT None

      integer(kind=iwp), intent(in):: NAS, iLo, iHi, jLo, jHi, LDS
      real(kind=wp), intent(In) :: SCA(NAS)
      real(kind=wp), intent(InOut)::  S(LDS,*)

      integer(kind=iwp) J, I

      DO J=jLo,jHi
        DO I=iLo,iHi
        S(I-iLo+1,J-jLo+1)=SCA(I)*SCA(J)*S(I-iLo+1,J-jLo+1)
        END DO
      END DO
      END SUBROUTINE S_SCALE

      SUBROUTINE V_SCALE (EIG,SCA,V,nRows,NAS,LDV,NIN,COND)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_module, only: ThrShS

      IMPLICIT None

      integer(kind=iwp), intent(in):: nRows, NAS, LDV, NIN
      real(kind=wp), intent(in):: EIG(NAS),SCA(NAS)
      real(kind=wp), Intent(inout):: V(LDV,*)
      real(kind=wp), Intent(out):: COND(NIN)

      integer(kind=iwp) jVEC, J, I, iVec
      real(kind=wp) EVAL, FACT, SZ

      jVEC=0
      DO J=1,NAS
        EVAL=EIG(J)
        IF(EVAL.GE.THRSHS) THEN
          jVEC=jVEC+1
          FACT=One/SQRT(EVAL)
          IF(jVEC.EQ.J) THEN
            CALL DSCAL_(nRows,FACT,V(1,J),1)
          ELSE
            CALL DYAX(nRows,FACT,V(1,J),1,V(1,jVEC),1)
          END IF
        END IF
      END DO
      IF (jVEC.NE.NIN) THEN
        WRITE(6,*) 'V_SCALE: '//
     &      'inconsitency in linear dependence removal, ABORT'
        call AbEnd()
      END IF
C Addition, for the scaled symmetric ON.
      DO I=1,nRows
        CALL DSCAL_(NIN,SCA(I),V(I,1),LDV)
      END DO
C The condition number, after scaling, disregarding linear dep.
      IF(NIN.GE.2) THEN
        DO jVEC=1,NIN
          SZ=Zero
          DO iVEC=1,nRows
            SZ=SZ+V(iVEC,jVEC)**2
          END DO
          COND(jVEC)=SZ
        END DO
      END IF
      END SUBROUTINE V_SCALE
#endif
