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
      SUBROUTINE SIGVEC_m(CIN,HC,HD,BM,SXN,G,H,DIA,F1,F2,X,C,NTRIAL)
C
C RASSCF Program version IBM-3090: SX section
C
C Purpose: Calculation of the SIGMA vector HC for a super-CI
C Hamiltonian.
C called from HMAT
CPAM01 Added miniscule constant times CIN to HC.
C
C ********** IBM-3090 Release 88 09 01 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CIN(*),HC(*),BM(*),SXN(*),G(*),H(*),DIA(*),F1(*),F2(*),
     &          X(*),C(*)
      DIMENSION HD(*)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='SIGVEC  ')
      Call qEnter('SIGVEC')
C Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
      NST=0
      DO 300 ITRIAL=1,NTRIAL
       NNST=NROOT+NST
C
C renormalize the C vector
C
       DO 10 I=1,NSXS
        C(I)=SXN(I)*CIN(I+NNST)
10     CONTINUE

C Remove any unwanted rotations from C:
       DO I=1,NSXS
        IF(HD(I+NROOT).GT.1.0D20) C(I)=0.0D0
       END DO

C Initialize sigma vector to zero.
       CALL DCOPY_(NROOT+NSXS,0.0D0,0,HC(NST+1),1)

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
C
C G-matrix contribution to HC (G*C)
C
        CALL DGEMM_('N','N',NIA,NAE,NIA,
     &              1.0D0,G(ISTIA),NIA,C(ISTBM+NST),NIA,
     &              1.0D0,HC(ISTBM+NNST),NIA)
C
C H-matrix contribution to HC (-C*H)
C
        IF(NAO.NE.0) CALL DGEMM_('N','N',NIA,NAO,NAE,
     &                           -1.0D0,C(ISTBM+NST),NIA,H(ISTH),NAE,
     &                           1.0D0,HC(ISTBM+NNST),NIA)
        IF(NAO*NEO.NE.0) CALL DGEMM_('N','T',NIA,NEO,NAO,
     &                         -1.0D0,C(ISTBM+NST),NIA,H(ISTH+NAO),NAE,
     &                          1.0D0,HC(ISTBM+NAO*NIA+NNST),NIA)
C
C First Fock matrix contribution D*C*FP
C
        CALL DGEMM_('N','N',
     &              NIA,NAE,NAE,
     &              1.0d0,C(ISTBM+NST),NIA,
     &              F2(ISTAE),NAE,
     &              0.0d0,X,NIA)
        CALL DGEMM_('N','N',NIA,NAE,NIA,
     &             1.0D0,DIA(ISTIA),NIA,X,NIA,
     &             1.0D0,HC(ISTBM+NNST),NIA)
C
C Second Fock matrix contribution FP*C*D
C
       IF(NAO.NE.0) THEN
         CALL DGEMM_('N','N',
     &               NIA,NAO,NIA,
     &               1.0d0,F1(ISTIA),NIA,
     &               C(ISTBM+NST),NIA,
     &               0.0d0,X,NIA)
         CALL DGEMM_('N','N',NIA,NAO,NAO,
     &              1.0D0,X,NIA,DIA(ISTIA+NIA*NIO+NIO),NIA,
     &              1.0D0,HC(ISTBM+NNST),NIA)
        ENDIF
C
98      ISTIA=ISTIA+NIA**2
        ISTAE=ISTAE+NAE**2
        ISTBM=ISTBM+NIA*NAE
        ISTH=ISTH+NAO*NAE
        ISTZ=ISTZ+(NAO**2-NAO)/2
100    CONTINUE
C
C Add diagonal contributions to the CI part
C
      IF(ICICP.NE.0) THEN
       DO 105 I=1,NROOT
       HC(NST+I)=HC(NST+I)+CIN(NST+I)*(ENER(I,ITER)-ENER(IROOT(1),ITER))
105    CONTINUE
      ENDIF
C
      If (NSXS.ne.0) Then
C
C BM contributions to CI part of sigma
C
*         CALL DGEMTX(NSXS,NROOT,1.0D0,BM,NSXS,C,1,HC(NST+1),1)
         CALL DGEMV_('T',NSXS,NROOT,1.0D0,BM,NSXS,C,1,1.0D0,HC(NST+1),1)
C
C BM contributions to SX part of sigma
C
         CALL DAXPY_(NSXS,CIN(NST+1),BM,1,HC(NNST+1),1)
C
      End If

C Remove any unwanted rotations:
       DO I=1,NSXS
        IF(HD(I+NROOT).GT.1.0D20) HC(I+NNST)=0.0D0
       END DO

C Adding a constant times C to the HC vectors at this point
C is equivalent to modifying the underlying minimization problem
C by adding a penalty for orbital rotations.
C This should be large enough that rotations get confined to values
C for which the Taylor expansion of rotations (which underlies the
C orbital optimization theory) is not too unrealistic.
C Example: if SXDAMP=0.0001, then a total rotation 'angle' of 0.1
C radians is penalized as an energy increase of 0.0001*(0.1)**2
C i.e., 10**(-6) a.u.
C Larger SXDAMP may make iterations sluggier.
C Too low SXDAMP may give erratic convergence in some exceptional
C cases.
C SXDAMP is set in READIN.
       CALL DAXPY_(NSXS,SXDAMP,C,1,HC(NNST+1),1)
C
C Renormalize the sigma vector
C
       DO 110 I=1,NSXS
        HC(I+NNST)=HC(I+NNST)*SXN(I)
110    CONTINUE
C
C Add Level shift part of the diagonal
C
      CALL DAXPY_(NSXS,SXSHFT,CIN(1+NNST),1,HC(1+NNST),1)
      NST=NST+NDIMSX
300   CONTINUE
C
C Test print out of the sigma vector
C
      IF(IPRLEV.GE.DEBUG) THEN
         Write(LF,1000) (HC(I),I=1,NDIMSX)
      END IF
1000  FORMAT(/1X,'Sigma vector in SIGVEC'/(1X,10F11.6))
C
      CALL QEXIT('SIGVEC')
      RETURN
      END
