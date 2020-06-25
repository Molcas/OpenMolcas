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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "para_info.fh"

      CALL QENTER('SBDIAG')

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
              CALL SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
            END IF
#else
            CALL SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
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

      CALL QEXIT('SBDIAG')

      RETURN
      END

      SUBROUTINE SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

* For fooling some compilers:
      DIMENSION WGRONK(2)

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

      SD = 0.0D0 ! dummy initialize

      CPU=0.0D0
      CONDNR=0.0D0
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS

      IF(NCOEF.EQ.0) RETURN

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(6,'("DEBUG> ",A12,A7,I2,A2,A6,A2,A5,I1)')
     &  'SBDIAG_SER: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM
      END IF

      NS=(NAS*(NAS+1))/2
      CALL GETMEM('LS','ALLO','REAL',LS,NS)
      IDS=IDSMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,2,WORK(LS),NS,IDS)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NS,WORK(LS),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ', FP
      END IF

C For some purposes, we need to save the diagonal elements:
      IF(BMATRIX.EQ.'YES') THEN
        IF(BTRANS.NE.'YES') THEN
          CALL GETMEM('LSD','ALLO','REAL',LSD,NAS)
          IDIAG=0
          DO I=1,NAS
            IDIAG=IDIAG+I
            WORK(LSD-1+I)=WORK(LS-1+IDIAG)
          END DO
        END IF
      END IF

C FIRST, FIND NIN ORTHONORMAL VECTORS BY SCALED SYMMETRIC ON.
C Addition, for the scaled symmetric ON: the S matrix is scaled
C to make the diagonal elements close to 1.
C Extremely small values give scale factor exactly zero.
      CALL GETMEM('LSCA','ALLO','REAL',LSCA,NAS)
      IDIAG=0
      DO I=1,NAS
        IDIAG=IDIAG+I
        SD=WORK(LS-1+IDIAG)
        IF(SD.GT.THRSHN) THEN
* Small variations of the scale factor were beneficial
          WORK(LSCA-1+I)=(1.0D00+DBLE(I)*3.0D-6)/SQRT(SD)
        ELSE
          WORK(LSCA-1+I)=0.0D0
        END IF
      END DO
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          WORK(LS-1+IJ)=WORK(LS-1+IJ)*
     &        WORK(LSCA-1+I)*WORK(LSCA-1+J)
        END DO
      END DO
C End of addition.
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NS,WORK(LS),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ', FP
      END IF

C DIAGONALIZE THE SCALED S MATRIX:
      CALL GETMEM('LVEC','ALLO','REAL',LVEC,NAS**2)
      CALL GETMEM('LEIG','ALLO','REAL',LEIG,NAS)

      CALL TIMING(CPU1,CPUE,TIO,TIOE)
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          WORK(LVEC-1+NAS*(J-1)+I)=WORK(LS-1+IJ)
        END DO
      END DO
      INFO=0
      call dsyev_('V','L',NAS,WORK(LVEC),NAS,WORK(LEIG),WGRONK,-1,INFO)
      NSCRATCH=INT(WGRONK(1))
      CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
      call dsyev_('V','U',NAS,WORK(LVEC),NAS,WORK(LEIG),WORK(LSCRATCH),
     &            NSCRATCH,INFO)
      CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
      CALL GETMEM('LS','FREE','REAL',LS,NS)

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1
      ! fingerprint eigenvalues
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS,WORK(LEIG),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT EIGENVALUE NORM: ', FP
      END IF

C Form orthonormal vectors by scaling eigenvectors
      NIN=0
      DO I=1,NAS
        EVAL=WORK(LEIG-1+I)
        IF(EVAL.LT.THRSHS) CYCLE
        FACT=1.0D00/SQRT(EVAL)
        NIN=NIN+1
        LVSTA=LVEC+NAS*(I-1)
        IF(NIN.EQ.I) THEN
          CALL DSCAL_(NAS,FACT,WORK(LVSTA),1)
        ELSE
          LVNEW=LVEC+NAS*(NIN-1)
          CALL DYAX(NAS,FACT,WORK(LVSTA),1,WORK(LVNEW),1)
        END IF
      END DO
      NINDEP(ISYM,ICASE)=NIN
      CALL GETMEM('LEIG','FREE','REAL',LEIG,NAS)
C Addition, for the scaled symmetric ON.
      DO I=1,NAS
        SCA=WORK(LSCA-1+I)
        CALL DSCAL_(NIN,SCA,WORK(LVEC-1+I),NAS)
      END DO
      CALL GETMEM('LSCA','FREE','REAL',LSCA,NAS)
C The condition number, after scaling, disregarding linear dep.
      IF(NIN.GE.2) THEN
        SZMIN=1.0D99
        SZMAX=0.0D0
        DO I=1,NIN
          SZ=DNRM2_(NAS,WORK(LVEC+NAS*(I-1)),1)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=(SZMAX/SZMIN)**2
      END IF
C End of addition.
      IF(NIN.EQ.0) THEN
        CALL GETMEM('LVEC','FREE','REAL',LVEC,NAS**2)
        RETURN
      END IF

      IF(BMATRIX.EQ.'NO') THEN
C In some calculations, we do not use B matrices.
C Just write the transformation matrix and branch out:
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LVEC),NAS*NIN,IDT)
        CALL GETMEM('LVEC','FREE','REAL',LVEC,NAS**2)
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
        CALL GETMEM('LBD','ALLO','REAL',LBD,NAS)
        NB=(NAS*(NAS+1))/2
        CALL GETMEM('LB','ALLO','REAL',LB,NB)
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,WORK(LB),NB,IDB)
        IDIAG=0
        DO I=1,NAS
          IDIAG=IDIAG+I
          WORK(LBD-1+I)=WORK(LB-1+IDIAG)
        END DO
        CALL GETMEM('LB','FREE','REAL',LB,NB)
C Now, the transformation matrix can be written out.
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LVEC),NAS*NIN,IDT)
        CALL GETMEM('LVEC','FREE','REAL',LVEC,NAS**2)
        DO I=1,NAS
          SD=WORK(LSD-1+I)+1.0d-15
          WORK(LBD-1+I)=WORK(LBD-1+I)/SD
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LBD),NAS,IDB)
        CALL GETMEM('LSD','FREE','REAL',LSD,NAS)
        CALL GETMEM('LBD','FREE','REAL',LBD,NAS)
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
      IDTMP=0
      CALL DDAFILE(LUSOLV,1,WORK(LVEC),NAS*NAUX,IDTMP)
      DO KSTA=NAUX+1,NIN,NAUX
        KEND=MIN(KSTA-1+NAUX,NIN)
        NCOL=1+KEND-KSTA
        LVSTA=LVEC+NAS*(KSTA-1)
        CALL DDAFILE(LUSOLV,1,WORK(LVSTA),NAS*NCOL,IDTMP)
      END DO
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS**2,WORK(LVEC),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)')
     &   'EIGENVECTOR NORM BEFORE B TRANS: ', FP
      END IF

C TRANSFORM B. NEW B WILL OVERWRITE AND DESTROY WORK(LVEC)
      IDB=IDBMAT(ISYM,ICASE)
      NB=NS
      CALL GETMEM('LB','ALLO','REAL',LB,NB)
      CALL DDAFILE(LUSBT,2,WORK(LB),NB,IDB)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NB,WORK(LB),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM: ', FP
      END IF

      CALL GETMEM('LBX','ALLO','REAL',LBX,NAS)
      CALL GETMEM('LXBX','ALLO','REAL',LXBX,NAS)
      DO J=NIN,1,-1
        LVSTA=LVEC+NAS*(J-1)
        CALL DCOPY_(NAS,[0.0D0],0,WORK(LBX),1)
#ifdef _CRAY_C90_
        CALL SSPMV('U',NAS,1.0D+00,WORK(LB),WORK(LVSTA),1,
     &                           1.0D+00,WORK(LBX),1)
#else
*        CALL DSLMX(NAS,1.0D+00,WORK(LB),WORK(LVSTA),1,
*     &                                   WORK(LBX),1)
        CALL DSPMV_('U',NAS,1.0D+00,WORK(LB),WORK(LVSTA),1,
     &                           1.0D+00,WORK(LBX),1)
#endif
C WORK(LBX): B * Vector number J.
        CALL DCOPY_(J,[0.0D0],0,WORK(LXBX),1)
        CALL DGEMM_('T','N',
     &              J,1,NAS,
     &              1.0d0,WORK(LVEC),NAS,
     &              WORK(LBX),NAS,
     &              0.0d0,WORK(LXBX),J)
C WORK(LXBX) CONTAINS NOW THE UPPERTRIANGULAR
C ELEMENTS OF THE J-th COLUMN OF TRANSFORMED B MATRIX.
        CALL DCOPY_(J,WORK(LXBX),1,WORK(LVSTA),1)
      END DO
      CALL GETMEM('LBX','FREE','REAL',LBX,NAS)
      CALL GETMEM('LXBX','FREE','REAL',LXBX,NAS)
      CALL GETMEM('LB','FREE','REAL',LB,NB)
C WORK(LVEC) HAS NOW BEEN DESTROYED (OVERWRITTEN BY NEW B).
C COPY TO TRIANGULAR STORAGE.
      NBNEW=(NIN*(NIN+1))/2
      CALL GETMEM('LB','ALLO','REAL',LB,NBNEW)
      DO J=1,NIN
        JOFF=(J*(J-1))/2
        CALL DCOPY_(J,WORK(LVEC+NAS*(J-1)),1,WORK(LB+JOFF),1)
      END DO
      CALL GETMEM('LVEC','FREE','REAL',LVEC,NAS**2)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NBNEW,WORK(LB),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

C DIAGONALIZE THE TRANSFORMED B MATRIX.
      CALL GETMEM('LEIG','ALLO','REAL',LEIG,NIN)
      CALL GETMEM('LVEC','ALLO','REAL',LVEC,NIN**2)
      CALL TIMING(CPU1,CPUE,TIO,TIOE)
C - Alt 0: Use diagonal approxim., if allowed:
      IF(BSPECT.NE.'YES')  THEN
        IDIAG=0
        DO I=1,NIN
          IDIAG=IDIAG+I
          WORK(LEIG-1+I)=WORK(LB-1+IDIAG)/SD
        END DO
      ELSE
        NBB=(NIN*(NIN+1))/2
        IJ=0
        DO J=1,NIN
          DO I=1,J
            IJ=IJ+1
            WORK(LVEC-1+NIN*(J-1)+I)=WORK(LB-1+IJ)
          END DO
        END DO
        CALL DSYEV_('V','U',NIN,WORK(LVEC),NIN,WORK(LEIG),
     &              WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
        CALL DSYEV_('V','U',NIN,WORK(LVEC),NIN,WORK(LEIG),
     &              WORK(LSCRATCH),NSCRATCH,INFO)
        CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
        CALL GETMEM('LB','FREE','REAL',LB,NB)
      END IF
      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NIN,WORK(LEIG),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

C The eigenvalues are written back at same position as the
C original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,WORK(LEIG),NIN,IDB)
      CALL GETMEM('LEIG','FREE','REAL',LEIG,NIN)

C Finally, we must form the composite transformation matrix,
C  = (NAS*NIN matrix on disk) * (NIN*NIN matrix in core).
C Assume enough space since we got rid of S/B matrices.
C Specifically, assume we have enough space for the two
C full matrices, plus an additional 19 columns of results.
      NAUX=MIN(19,NIN)
      CALL GETMEM('LTRANS','ALLO','REAL',LTRANS,NAS*NIN)
      CALL GETMEM('LAUX'  ,'ALLO','REAL',LAUX  ,NAS*NAUX)
      IDTMP=0
      CALL DDAFILE(LUSOLV,2,WORK(LAUX),NAS*NAUX,IDTMP)
      IF(BTRANS.EQ.'YES') THEN
        CALL DGEMM_('N','N',
     &              NAS,NIN,NAUX,
     &              1.0d0,WORK(LAUX),NAS,
     &              WORK(LVEC),NIN,
     &              0.0d0,WORK(LTRANS),NAS)
      ELSE
        CALL DCOPY_(NAS*NAUX,WORK(LAUX),1,WORK(LTRANS),1)
      END IF
      DO KSTA=NAUX+1,NIN,NAUX
        KEND=MIN(KSTA-1+NAUX,NIN)
        NCOL=1+KEND-KSTA
        CALL DDAFILE(LUSOLV,2,WORK(LAUX),NAS*NCOL,IDTMP)
        IF(BTRANS.EQ.'YES') THEN
          CALL DGEMM_('N','N',NAS,NIN,NCOL,1.0D00,
     &              WORK(LAUX),NAS,WORK(LVEC-1+KSTA),NIN,
     &              1.0D00,WORK(LTRANS),NAS)
        ELSE
          LTRANS1=LTRANS+NAS*(KSTA-1)
          CALL DCOPY_(NAS*NCOL,WORK(LAUX),1,WORK(LTRANS1),1)
        END IF
      END DO
      CALL GETMEM('LAUX'  ,'FREE','REAL',LAUX  ,NAS*NAUX)
      CALL GETMEM('LVEC','FREE','REAL',LVEC,NIN**2)
      IDT=IDTMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,WORK(LTRANS),NAS*NIN,IDT)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS*NIN,WORK(LTRANS),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', FP
      END IF

C-SVC: compute S*T and store on disk for later use by RHS vector
C      utilities.
      NS=(NAS*(NAS+1))/2
      CALL GETMEM('LS','ALLO','REAL',LS,NS)
      IDS=IDSMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,2,WORK(LS),NS,IDS)
      CALL GETMEM('LST','ALLO','REAL',LST,NAS*NIN)
      CALL DCOPY_(NAS*NIN,[0.0D0],0,WORK(LST),1)
      CALL TRIMUL(NAS,NIN,1.0D00,WORK(LS),WORK(LTRANS),
     &            NAS,WORK(LST),NAS)
      CALL GETMEM('LS','FREE','REAL',LS,NS)
      CALL GETMEM('LTRANS','FREE','REAL',LTRANS,NAS*NIN)
      IDST=IDSTMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,WORK(LST),NAS*NIN,IDST)
      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS*NIN,WORK(LST),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'STMAT NORM: ', FP
      END IF
      CALL GETMEM('LST','FREE','REAL',LST,NAS*NIN)

      END

C SVC2010: The following subroutine computes the transformation
C matrices using global arrays rather than replicate arrays.  Currently,
C the maximum memory scales approximately as 3*(NAS**2), which is twice
C what is needed in the replicate routines.  It can be reduced to
C 2*(NAS**2), if the matrix multiplications are rewritten to occur in
C batch mode.  However, unlike in the replicate routine, this amount is
C divided over processors.
#ifdef _MOLCAS_MPP_
      SUBROUTINE SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

C-SVC20100902: global arrays header files
#include "global.fh"
#include "mafdecls.fh"
#ifndef SCALAPACK
      DIMENSION WGRONK(2)
#endif
      LOGICAL bSTAT
      CHARACTER(2) cSYM,cCASE
      LOGICAL KING

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

      CPU=0.0D0
      CONDNR=0.0D0
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
        CALL GETMEM('COL','ALLO','REAL',LCOL,NCOL)
        CALL GETMEM('TMP','ALLO','REAL',LTMP,NTMP)
        iOFF=0
        DO J=1,NAS
          call GA_Get (lg_S, 1, J, J, J, WORK(LCOL), NAS)
          CALL DCOPY_(J,WORK(LCOL),1,WORK(LTMP+iOFF),1)
          iOFF=iOFF+J
        END DO
        CALL GETMEM('COL','FREE','REAL',LCOL,NCOL)
        IDS=IDSMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LTMP),NTMP,IDS)
        CALL GETMEM('TMP','FREE','REAL',LTMP,NTMP)
      END IF

C Save the diagonal elements from the S matrix for easy access later on.
C FIXME: nicer way to do this?
      CALL GETMEM('SD','ALLO','REAL',LSD,NAS)
      CALL DCOPY_(NAS,[0.0D0],0,WORK(LSD),1)
      myRank = GA_NodeID()
      call GA_Distribution (lg_S, myRank, iLo, iHi, jLo, jHi)
      ISTA=MAX(ILO,JLO)
      IEND=MIN(IHI,JHI)
      IF (ISTA.NE.0) THEN
        call GA_Access (lg_S, iLo, iHi, jLo, jHi, mS, LDS)
        DO I=ISTA,IEND
          WORK(LSD+I-1)=DBL_MB(mS+I-ILO+LDS*(I-JLO))
        END DO
        call GA_Release (lg_S, iLo, iHi, jLo, jHi)
      END IF
      CALL GADSUM (WORK(LSD),NAS)

C Calculate the scaling factors and store them in array LSCA.
      CALL GETMEM('SCA','ALLO','REAL',LSCA,NAS)
      DO I=1,NAS
        SD=WORK(LSD+I-1)
        IF(SD.GT.THRSHN) THEN
          WORK(LSCA-1+I)=(1.0D00+DBLE(I)*3.0D-6)/SQRT(SD)
        ELSE
          WORK(LSCA-1+I)=0.0D0
        END IF
      END DO

C Scale the elements S(I,J) with the factor SCA(I)*SCA(J).
      myRank = GA_NodeID()
      call GA_Distribution (lg_S, myRank, iLo, iHi, jLo, jHi)
      IF (iLo.NE.0) THEN
        call GA_Access (lg_S, iLo, iHi, jLo, jHi, mS, LDS)
        call S_SCALE (NAS,WORK(LSCA),DBL_MB(mS),iLo,iHi,jLo,jHi,LDS)
        call GA_Release_Update (lg_S, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()
      CALL TIMING(CPU1,CPUE,TIO,TIOE)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_S,NAS)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ', FP
      END IF

      CALL GETMEM('LEIG','ALLO','REAL',LEIG,NAS)
      CALL DCOPY_(NAS,[0.0D0],0,WORK(LEIG),1)
C Diagonalize the global array S.  Some (old) reports about parallel
C performance recommend PDSYEVX or PDSYEVR as fastest methods if
C eigenvectors are needed (FIXME: should time this).  For the linear
C dependence removal, split eigenvectors in horizontal stripes so that
C each processor has a row window of all column vectors
#ifdef SCALAPACK
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      CALL GA_PDSYEVX_ (lg_S, lg_V, WORK(LEIG), 0)
      bSTAT = GA_Destroy (lg_S)
#else
C here for the non-ScaLAPACK version: copy matrix to master process,
C diagonalize using the serial DSYEV routine, and copy the resulting
C eigenvectors back to a global array.  Then distribute the eigenvalues.
      IF (myRank.EQ.0) THEN
        CALL GETMEM('LVEC','ALLO','REAL',LVEC,NAS**2)
        CALL GA_Get (lg_S, 1, NAS, 1, NAS, WORK(LVEC), NAS)
      END IF
      bSTAT = GA_Destroy (lg_S)
      IF (myRank.EQ.0) THEN
        CALL DSYEV_('V','L',NAS,WORK(LVEC),NAS,WORK(LEIG),
     &              WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
        CALL DSYEV_('V','L',NAS,WORK(LVEC),NAS,WORK(LEIG),
     &              WORK(LSCRATCH),NSCRATCH,INFO)
        CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
      END IF
      CALL PSBMAT_GETMEM('VMAT',lg_V,NAS)
      IF (myRank.EQ.0) THEN
        CALL GA_Put (lg_V, 1, NAS, 1, NAS, WORK(LVEC), NAS)
        CALL GETMEM('LVEC','FREE','REAL',LVEC,NAS**2)
      END IF
      CALL GADSUM(WORK(LEIG),NAS)
#endif

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NAS,WORK(LEIG),1)
        WRITE(6,'("DEBUG> ",A,ES21.14)') 'SMAT EIGENVALUE NORM: ', FP
      END IF

      NIN=0
      DO J=1,NAS
        IF(WORK(LEIG+J-1).GE.THRSHS) NIN=NIN+1
      END DO
      NINDEP(ISYM,ICASE)=NIN
      IF (NIN.EQ.0) THEN
        CALL GETMEM('SCA','FREE','REAL',LSCA,NAS)
        CALL GETMEM('LEIG','FREE','REAL',LEIG,NAS)
        CALL GETMEM('SD','FREE','REAL',LSD,NAS)
        bSTAT = GA_Destroy (lg_V)
        RETURN
      END IF

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1

      CALL GETMEM('COND','ALLO','REAL',LCOND,NIN)
      CALL DCOPY_(NIN,[0.0D0],0,WORK(LCOND),1)
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
        call V_SCALE (WORK(LEIG),WORK(LSCA+iLo-1),DBL_MB(mV),
     &              iHi-iLo+1,jHi-jLo+1,LDV,NIN,WORK(LCOND))
        call GA_Release_Update (lg_V, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()

      CALL GETMEM('LEIG','FREE','REAL',LEIG,NAS)
      CALL GETMEM('LSCA','FREE','REAL',LSCA,NAS)

C The condition number, after scaling, disregarding linear dep.
C FIXME: adapt to local subroutine for global array lg_V
      IF(NIN.GE.2) THEN
        CALL GADSUM (WORK(LCOND),NIN)
        SZMIN=1.0D99
        SZMAX=0.0D0
        DO I=1,NIN
          SZ=WORK(LCOND+I-1)
          SZMIN=MIN(SZMIN,SZ)
          SZMAX=MAX(SZMAX,SZ)
        END DO
        CONDNR=SZMAX/SZMIN
      END IF
      CALL GETMEM('COND','FREE','REAL',LCOND,NIN)

C Copy the NIN non-linear dependent eigenvectors to the transformation
C matrix T(NAS,NIN).
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
#ifdef _GA_
      call GA_Copy_Patch ('N', lg_V, 1, NAS, 1, NIN,
     &                         lg_T, 1, NAS, 1, NIN)
#else
      call GA_Distribution (lg_V, myRank, iLoV, iHiV, jLoV, jHiV)
      call GA_Distribution (lg_T, myRank, iLoT, iHiT, jLoT, jHiT)
      IF (iLoV.NE.0 .AND. iLoT.NE.0) THEN
        NROW=iHiT-iLoT+1
        NCOL=jHiT-jLoT+1
        call GA_Access (lg_V, iLoV, iHiV, jLoV, jHiV, mV, LDV)
        call GA_Access (lg_T, iLoT, iHiT, jLoT, jHiT, mT, LDT)
        IF (LDV.NE.LDT) THEN
          WRITE(6,'(1X,A)') 'SBDIAG_MPP: LDV != LDT, abort!'
        END IF
        CALL DCOPY_(NROW*NCOL,DBL_MB(mV),1,DBL_MB(mT),1)
        call GA_Release_Update (lg_T, iLoT, iHiT, jLoT, jHiT)
        call GA_Release (lg_V, iLoV, iHiV, jLoV, jHiV)
      END IF
#endif
      bStat = GA_Destroy (lg_V)

      IF(BMATRIX.EQ.'NO      ') THEN
C In some calculations, we do not use B matrices.
C Write the T matrix to disk and exit.  FIXME: This
C should be removed when the transformation matrices are stored as disk
C resident arrays only.
        IF (KING()) THEN
          CALL GETMEM('LTRANS','ALLO','REAL',LTRANS,NAS*NIN)
          call GA_Get (lg_T, 1, NAS, 1, NIN, WORK(LTRANS), NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,WORK(LTRANS),1)
            WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,WORK(LTRANS),NAS*NIN,IDT)
          CALL GETMEM('LTRANS','FREE','REAL',LTRANS,NAS*NIN)
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
        CALL GETMEM('BD','ALLO','REAL',LBD,NAS)
        CALL DCOPY_(NAS,[0.0D0],0,WORK(LBD),1)
        myRank = GA_NodeID()
        call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
        ISTA=MAX(ILO,JLO)
        IEND=MIN(IHI,JHI)
        IF (ISTA.NE.0) THEN
          call GA_Access (lg_B, iLo, iHi, jLo, jHi, mB, LDB)
          DO I=ISTA,IEND
            WORK(LBD+I-1)=DBL_MB(mB+I-ILO+LDB*(I-JLO))
          END DO
          call GA_Release (lg_B, iLo, iHi, jLo, jHi)
        END IF
        CALL GADSUM (WORK(LBD),NAS)
        bStat = GA_Destroy (lg_B)
        DO I=1,NAS
          SD=WORK(LSD-1+I)+1.0d-15
          WORK(LBD-1+I)=WORK(LBD-1+I)/SD
        END DO
        IDB=IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LBD),NAS,IDB)
        CALL GETMEM('BD','FREE','REAL',LBD,NAS)
C Write the transformation matrices after the diagonal values of B.
C FIXME: This should be removed when the transformation matrices are
C stored as disk resident arrays only.
        IF (KING()) THEN
          CALL GETMEM('LTRANS','ALLO','REAL',LTRANS,NAS*NIN)
          call GA_Get (lg_T, 1, NAS, 1, NIN, WORK(LTRANS), NAS)
          IF (iPrGlb.GE.INSANE) THEN
            dTRANS=dNRM2_(NAS*NIN,WORK(LTRANS),1)
            WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
          END IF
          IDT=IDTMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,WORK(LTRANS),NAS*NIN,IDT)
          CALL GETMEM('LTRANS','FREE','REAL',LTRANS,NAS*NIN)
        END IF
        bStat = GA_Destroy (lg_T)
        RETURN
      END IF
      CALL GETMEM('SD','FREE','REAL',LSD,NAS)

C TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_T,NAS)
        WRITE(6,'(1X,A,ES21.14)')
     &   'EIGENVECTOR NORM BEFORE B TRANS: ', FP
      END IF

      CALL PSBMAT_GETMEM ('B',lg_B,NAS)
      CALL PSBMAT_READ ('B',iCase,iSym,lg_B,NAS)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_B,NAS)
        WRITE(6,'(1X,A,ES21.14)') 'BMAT NORM: ', FP
      END IF

C FIXME: Perform transformation of B using horizontal stripes of B or
C vertical stripes of T to reduce memory usage if necessary as indicated
C by the available memory, which is now scaling as approx. 3*(NAS**2).
#ifdef _GA_
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, 1.0D0,
     &               lg_B, lg_T, 0.0D0, lg_X )
      bStat = GA_Destroy (lg_B)

      CALL GA_CREATE_STRIPED ('H',NIN,NIN,'BMAT',lg_B)
      call GA_DGEMM ('T', 'N', NIN, NIN, NAS, 1.0D0,
     &               lg_T, lg_X, 0.0D0, lg_B )
      bStat = GA_Destroy (lg_X)
#else
C FIXME: Remove second part if DGA supports DGEMM.
      IF (KING()) THEN
        CALL GETMEM('TMAT','ALLO','REAL',LT,NAS*NIN)
        CALL GA_GET (LG_T, 1, NAS, 1, NIN, WORK(LT), NAS)
        CALL GETMEM('BMAT','ALLO','REAL',LB,NAS*NAS)
        CALL GA_Get (lg_B, 1, NAS, 1, NAS, WORK(LB), NAS)
        CALL GETMEM('XMAT','ALLO','REAL',LX,NAS*NIN)
        CALL DGEMM_ ('N', 'N', NAS, NIN, NAS, 1.0D0,
     &               WORK(LB), NAS, WORK(LT), NAS,
     &               0.0D0, WORK(LX), NAS)
        CALL GETMEM('BMAT','FREE','REAL',LB,NAS*NAS)
      END IF
      bStat = GA_Destroy (lg_B)
      CALL PSBMAT_GETMEM ('B',lg_B,NIN)
      IF (KING()) THEN
        CALL GETMEM('BMAT','ALLO','REAL',LB,NIN*NIN)
        CALL DGEMM_ ('T', 'N', NIN, NIN, NAS, 1.0D0,
     &               WORK(LT), NAS, WORK(LX), NAS,
     &               0.0D0, WORK(LB), NIN)
        CALL GA_Put (lg_B, 1, NIN, 1, NIN, WORK(LB), NIN)
        CALL GETMEM('BMAT','FREE','REAL',LB,NIN*NIN)
        CALL GETMEM('XMAT','FREE','REAL',LX,NAS*NIN)
        CALL GETMEM('TMAT','FREE','REAL',LT,NAS*NIN)
      END IF
#endif
      bStat = GA_Destroy (lg_T)

      IF (IPRGLB.GE.INSANE) THEN
        FP=PSBMAT_FPRINT(lg_B,NIN)
        WRITE(6,'(1X,A,ES21.14)') 'BMAT NORM AFTER TRANS: ', FP
      END IF

      CALL TIMING(CPU1,CPUE,TIO,TIOE)

C Diagonalize the transformed B matrix.
      CALL GETMEM('LEIG','ALLO','REAL',LEIG,NIN)
      CALL DCOPY_(NIN,[0.0D0],0,WORK(LEIG),1)
      IF(BSPECT.NE.'YES')  THEN
C Use diagonal approxim., if allowed.
C        call GA_Fill (lg_V, 0.0D0)
        call GA_Zero (lg_V)
C FIXME: this original code seemed wrong, using uninitialized SD?
*       IDIAG=1
*       DO I=1,NIN
*         WORK(LEIG-1+I)=WORK(LB-1+IDIAG)/SD
*         IDIAG=IDIAG+1+NIN-I
*       END DO
        WRITE(6,*) 'GLOB_SBDIAG: option not implemented'
        call AbEnd()
      ELSE
#ifdef SCALAPACK
        CALL GA_CREATE_STRIPED ('H',NIN,NIN,'VMAT',lg_V)
        CALL GA_PDSYEVX_ (lg_B, lg_V, WORK(LEIG), 0)
        bStat = GA_Destroy (lg_B)
#else
        IF (myRank.EQ.0) THEN
          CALL GETMEM('LVEC','ALLO','REAL',LVEC,NIN**2)
          CALL GA_Get (lg_B, 1, NIN, 1, NIN, WORK(LVEC), NIN)
        END IF
        bSTAT = GA_Destroy (lg_B)
        IF (myRank.EQ.0) THEN
          call dsyev_('V','L',NIN,WORK(LVEC),NIN,WORK(LEIG),
     &               WGRONK,-1,INFO)
          NSCRATCH=INT(WGRONK(1))
          CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
          call dsyev_('V','L',NIN,WORK(LVEC),NIN,WORK(LEIG),
     &               WORK(LSCRATCH),NSCRATCH,INFO)
          CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
        END IF
        call GA_Sync()
        CALL GA_CREATE_STRIPED ('H',NIN,NIN,'VMAT',lg_V)
        IF (myRank.EQ.0) THEN
          CALL GA_Put (lg_V, 1, NIN, 1, NIN, WORK(LVEC), NIN)
          CALL GETMEM('LVEC','FREE','REAL',LVEC,NIN**2)
        END IF
        CALL GADSUM(WORK(LEIG),NIN)
#endif
      END IF

      IF (IPRGLB.GE.INSANE) THEN
        FP=DNRM2_(NIN,WORK(LEIG),1)
        WRITE(6,'(1X,A,ES21.14)') 'BMAT EIGENVALUE NORM: ', FP
      END IF

C The eigenvalues are written back at same position as the
C original B matrix, which is destroyed:
      IDB=IDBMAT(ISYM,ICASE)
      CALL DDAFILE(LUSBT,1,WORK(LEIG),NIN,IDB)
      CALL GETMEM('LEIG','FREE','REAL',LEIG,NIN)

      CALL TIMING(CPU2,CPUE,TIO,TIOE)
      CPU=CPU+CPU2-CPU1

C Finally, we must form the composite transformation matrix: T(NAS,NIN)
C matrix on disk * V(NIN*NIN) matrix in core.  FIXME: for now, asume
C there is enough memory for the full transformation, scaling as
C approx. 3*(NAS**2).  Should be determined by the available memory.
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'XMAT',lg_X)
      CALL PSBMAT_READ ('T',iCase,iSym,lg_X,NAS*NIN)
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
#ifdef _GA_
      call GA_DGEMM ('N', 'N', NAS, NIN, NIN, 1.0D0,
     &               lg_X, lg_V, 0.0D0, lg_T )
#else
      IF (KING()) THEN
        CALL GETMEM('XMAT','ALLO','REAL',LX,NAS*NIN)
        CALL GA_GET (LG_X, 1, NAS, 1, NIN, WORK(LX), NAS)
        CALL GETMEM('VMAT','ALLO','REAL',LV,NIN*NIN)
        CALL GA_Get (lg_V, 1, NIN, 1, NIN, WORK(LV), NIN)
        CALL GETMEM('TMAT','ALLO','REAL',LT,NAS*NIN)
        CALL DGEMM_ ('N', 'N', NAS, NIN, NIN, 1.0D0,
     &               WORK(LX), NAS, WORK(LV), NIN,
     &               0.0D0, WORK(LT), NAS)
        CALL GA_Put (lg_T, 1, NAS, 1, NIN, WORK(LT), NAS)
        CALL GETMEM('XMAT','FREE','REAL',LX,NAS*NIN)
        CALL GETMEM('VMAT','FREE','REAL',LV,NIN*NIN)
      END IF
#endif
      bStat = GA_Destroy (lg_X)
      bStat = GA_Destroy (lg_V)

C Write the composite transformation matrix to disk.
      CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)

C Additonally, compute S*T and store in on disk for later use by the RHS
C vector utitlities
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)

      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'STMAT',lg_ST)
#ifdef _GA_
      call GA_DGEMM ('N', 'N', NAS, NIN, NAS, 1.0D0,
     &               lg_S, lg_T, 0.0D0, lg_ST )
#else
      IF (KING()) THEN
        CALL GETMEM('SMAT','ALLO','REAL',LS,NAS*NAS)
        CALL GA_GET (LG_S, 1, NAS, 1, NAS, WORK(LS), NAS)
        CALL GETMEM('STMAT','ALLO','REAL',LST,NAS*NIN)
        CALL DGEMM_ ('N', 'N', NAS, NIN, NAS, 1.0D0,
     &               WORK(LS), NAS, WORK(LT), NAS,
     &               0.0D0, WORK(LST), NAS)
        CALL GA_Put (lg_ST, 1, NAS, 1, NIN, WORK(LST), NAS)
        CALL GETMEM('STMAT','FREE','REAL',LST,NAS*NIN)
        CALL GETMEM('SMAT','FREE','REAL',LS,NAS*NAS)
        CALL GETMEM('TMAT','FREE','REAL',LT,NAS*NIN)
      END IF
#endif
      bStat = GA_Destroy (lg_S)

      CALL PSBMAT_WRITE ('M',iCase,iSym,lg_ST,NAS*NIN)
      bStat = GA_Destroy (lg_ST)

C For now, also keep the transformation matrix on disk as a
C replicate array.  FIXME: Should be removed later.
      IF (KING()) THEN
        CALL GETMEM('LTRANS','ALLO','REAL',LTRANS,NAS*NIN)
        call GA_Get (lg_T, 1, NAS, 1, NIN, WORK(LTRANS), NAS)
        dTRANS=dNRM2_(NAS*NIN,WORK(LTRANS),1)
        IF (iPrGlb.GE.INSANE) THEN
          WRITE(6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ', dTRANS
        END IF
        IDT=IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WORK(LTRANS),NAS*NIN,IDT)
        CALL GETMEM('LTRANS','FREE','REAL',LTRANS,NAS*NIN)
      END IF

      call ga_sync()
      bStat = GA_Destroy (lg_T)

      END

      SUBROUTINE S_SCALE (NAS,SCA,S,iLo,iHi,jLo,jHi,LDS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      DIMENSION SCA(NAS),S(LDS,*)
      DO J=jLo,jHi
        DO I=iLo,iHi
        S(I-iLo+1,J-jLo+1)=SCA(I)*SCA(J)*S(I-iLo+1,J-jLo+1)
        END DO
      END DO
      END

      SUBROUTINE V_SCALE (EIG,SCA,V,nRows,NAS,LDV,NIN,COND)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      DIMENSION EIG(NAS),SCA(NAS),V(LDV,*),COND(NIN)
      jVEC=0
      DO J=1,NAS
        EVAL=EIG(J)
        IF(EVAL.GE.THRSHS) THEN
          jVEC=jVEC+1
          FACT=1.0D00/SQRT(EVAL)
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
          SZ=0.0D0
          DO iVEC=1,nRows
            SZ=SZ+V(iVEC,jVEC)**2
          END DO
          COND(jVEC)=SZ
        END DO
      END IF
      END
#endif
