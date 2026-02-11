!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE SIGVEC(CIN,HC,HD,BM,SXN,G,H,DIA,F1,F2,X,C,NTRIAL)
!
! RASSCF Program version IBM-3090: SX section
!
! Purpose: Calculation of the SIGMA vector HC for a super-CI
! Hamiltonian.
! called from HMAT
!PAM01 Added miniscule constant times CIN to HC.
!
! ********** IBM-3090 Release 88 09 01 **********
!
      use rasscf_global, only: ICICP, ITER, NDIMSX, NROOT, NSXS,        &
     &                         SXSHFT, IROOT, ENER
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NASH,NISH,NSSH,SXDAMP

      IMPLICIT NONE
      REAL*8 CIN(*),HC(*),BM(*),SXN(*),G(*),H(*),DIA(*),F1(*),F2(*),    &
     &          X(*),C(*)
      INTEGER NTRIAL

      Character(LEN=16), Parameter :: ROUTINE='SIGVEC  '
      REAL*8 HD(*)
      INTEGER :: I, iPrLev, ISTAE, ISTBM, ISTH, ISTIA, ISTZ, ISYM,      &
     &           ITRIAL, NAE, NAO, NEO, NIA, NIO, NNST, NST

! Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
      NST=0
      DO 300 ITRIAL=1,NTRIAL
       NNST=NROOT+NST
!
! renormalize the C vector
!
       DO 10 I=1,NSXS
        C(I)=SXN(I)*CIN(I+NNST)
10     CONTINUE

! Remove any unwanted rotations from C:
       DO I=1,NSXS
        IF(HD(I+NROOT).GT.1.0D20) C(I)=0.0D0
       END DO

! Initialize sigma vector to zero.
       CALL DCOPY_(NROOT+NSXS,[0.0D0],0,HC(NST+1),1)

       ISTIA=1
       ISTAE=1
       ISTBM=1
       ISTH=1
       ISTZ=0
       DO 100 ISYM=1,NSYM
        NIO=NISH(ISYM)
        NAO=NASH(ISYM)
        NEO=NSSH(ISYM)
        NIA=NIO+NAO
        NAE=NAO+NEO
        IF(NIA.EQ.0.OR.NAE.EQ.0) GO TO 98
!
! G-matrix contribution to HC (G*C)
!
        CALL DGEMM_('N','N',NIA,NAE,NIA,                                &
     &              1.0D0,G(ISTIA),NIA,C(ISTBM+NST),NIA,                &
     &              1.0D0,HC(ISTBM+NNST),NIA)
!
! H-matrix contribution to HC (-C*H)
!
        IF(NAO.NE.0) CALL DGEMM_('N','N',NIA,NAO,NAE,                   &
     &                           -1.0D0,C(ISTBM+NST),NIA,H(ISTH),NAE,   &
     &                           1.0D0,HC(ISTBM+NNST),NIA)
        IF(NAO*NEO.NE.0) CALL DGEMM_('N','T',NIA,NEO,NAO,               &
     &                         -1.0D0,C(ISTBM+NST),NIA,H(ISTH+NAO),NAE, &
     &                          1.0D0,HC(ISTBM+NAO*NIA+NNST),NIA)
!
! First Fock matrix contribution D*C*FP
!
        CALL DGEMM_('N','N',                                            &
     &              NIA,NAE,NAE,                                        &
     &              1.0d0,C(ISTBM+NST),NIA,                             &
     &              F2(ISTAE),NAE,                                      &
     &              0.0d0,X,NIA)
        CALL DGEMM_('N','N',NIA,NAE,NIA,                                &
     &             1.0D0,DIA(ISTIA),NIA,X,NIA,                          &
     &             1.0D0,HC(ISTBM+NNST),NIA)
!
! Second Fock matrix contribution FP*C*D
!
       IF(NAO.NE.0) THEN
         CALL DGEMM_('N','N',                                           &
     &               NIA,NAO,NIA,                                       &
     &               1.0d0,F1(ISTIA),NIA,                               &
     &               C(ISTBM+NST),NIA,                                  &
     &               0.0d0,X,NIA)
         CALL DGEMM_('N','N',NIA,NAO,NAO,                               &
     &              1.0D0,X,NIA,DIA(ISTIA+NIA*NIO+NIO),NIA,             &
     &              1.0D0,HC(ISTBM+NNST),NIA)
        ENDIF
!
98      ISTIA=ISTIA+NIA**2
        ISTAE=ISTAE+NAE**2
        ISTBM=ISTBM+NIA*NAE
        ISTH=ISTH+NAO*NAE
        ISTZ=ISTZ+(NAO**2-NAO)/2
100    CONTINUE
!
! Add diagonal contributions to the CI part
!
      IF(ICICP.NE.0) THEN
       DO 105 I=1,NROOT
       HC(NST+I)=HC(NST+I)+CIN(NST+I)*(ENER(I,ITER)-ENER(IROOT(1),ITER))
105    CONTINUE
      ENDIF
!
      If (NSXS.ne.0) Then
!
! BM contributions to CI part of sigma
!
!         CALL DGEMTX(NSXS,NROOT,1.0D0,BM,NSXS,C,1,HC(NST+1),1)
         CALL DGEMV_('T',NSXS,NROOT,1.0D0,BM,NSXS,C,1,1.0D0,HC(NST+1),1)
!
! BM contributions to SX part of sigma
!
         CALL DAXPY_(NSXS,CIN(NST+1),BM,1,HC(NNST+1),1)
!
      End If

! Remove any unwanted rotations:
       DO I=1,NSXS
        IF(HD(I+NROOT).GT.1.0D20) HC(I+NNST)=0.0D0
       END DO

! Adding a constant times C to the HC vectors at this point
! is equivalent to modifying the underlying minimization problem
! by adding a penalty for orbital rotations.
! This should be large enough that rotations get confined to values
! for which the Taylor expansion of rotations (which underlies the
! orbital optimization theory) is not too unrealistic.
! Example: if SXDAMP=0.0001, then a total rotation 'angle' of 0.1
! radians is penalized as an energy increase of 0.0001*(0.1)**2
! i.e., 10**(-6) a.u.
! Larger SXDAMP may make iterations sluggier.
! Too low SXDAMP may give erratic convergence in some exceptional
! cases.
! SXDAMP is set in READIN.
       CALL DAXPY_(NSXS,SXDAMP,C,1,HC(NNST+1),1)
!
! Renormalize the sigma vector
!
       DO 110 I=1,NSXS
        HC(I+NNST)=HC(I+NNST)*SXN(I)
110    CONTINUE
!
! Add Level shift part of the diagonal
!
      CALL DAXPY_(NSXS,SXSHFT,CIN(1+NNST),1,HC(1+NNST),1)
      NST=NST+NDIMSX
300   CONTINUE
!
! Test print out of the sigma vector
!
      IF(IPRLEV.GE.DEBUG) THEN
         Write(LF,1000) (HC(I),I=1,NDIMSX)
      END IF
1000  FORMAT(/1X,'Sigma vector in SIGVEC'/(1X,10F11.6))
!
      RETURN
      END
