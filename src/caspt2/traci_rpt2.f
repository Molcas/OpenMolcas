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
      SUBROUTINE TRACI_RPT2(ISTART,NDIM,XMAT,STSYM,NCI,CI)
      use definitions, only: iwp, wp
      use constants, only: Zero, Half, One, OneHalf
      use gugx, only: LEVEL, SGS, CIS, EXS
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT None
      integer(kind=iwp), intent(in):: ISTART, NDIM, STSym, NCI
      real(kind=wp), intent(inout):: XMAT(NDIM,NDIM)
      real(kind=wp), intent(inout):: CI(NCI)

      real(kind=wp), ALLOCATABLE:: XSAV(:,:), TVEC(:), SGM(:)
      real(kind=wp) Fact, SCL, XJM
      integer(kind=iwp) I, IORB, J, JORB, LI, LJ, M

      IF (NDIM.LE.0) RETURN

      CALL mma_allocate(XSAV,NDIM,NDIM,Label='XSAV')
      XSAV(:,:)=XMAT(:,:)
      CALL mma_allocate(TVEC,NDIM,LABEL='TVEC')
      CALL mma_allocate(SGM,NCI,LABEL='SGM')
      SGM(:)=Zero

      DO J=1,NDIM
        FACT=One/XMAT(J,J)
        DO I=1,NDIM
          TVEC(I)=-FACT*XMAT(I,J)
          XMAT(I,J)=Zero
        END DO
        TVEC(J)=FACT
        XMAT(J,J)=One
C Array T now contains a factor of XMAT of the form
C (e(1),..e(k-1),T,..,e(n)), where e(i) is the standard
C unit column vector, with elements Kronecker(l,i).
C Apply its inverse to XMAT.
        DO M=J+1,NDIM
          XJM=XMAT(J,M)
          DO I=1,NDIM
            XMAT(I,M)=XMAT(I,M)+TVEC(I)*XJM
          END DO
          XMAT(J,M)=TVEC(J)*XJM
        END DO
C Transform CI array:
C CI:=( 1 + Sum(U(I)E(IJ)) + (1/2)Sum(U(I)U(M)E(IJ,MJ)) ) CI,
C where U(I) = T(I)-Kronecker(I,J).
        JORB=ISTART-1+J
        LJ=LEVEL(JORB)
        CALL DYAX(NCI,(OneHalf-Half*TVEC(J)),CI,1,SGM,1)
        DO I=1,NDIM
          IORB=ISTART-1+I
          LI=LEVEL(IORB)
          SCL=Half*TVEC(I)
          IF(I.EQ.J) SCL=SCL-Half
          CALL SIGMA1(SGS,CIS,EXS,LI,LJ,SCL,STSYM,CI,SGM)
        END DO
        DO I=1,NDIM
          IORB=ISTART-1+I
          LI=LEVEL(IORB)
          SCL=TVEC(I)
          IF(I.EQ.J) SCL=SCL-One
          CALL SIGMA1(SGS,CIS,EXS,LI,LJ,SCL,STSYM,SGM,CI)
        END DO

      End Do

      CALL mma_deallocate(SGM)
      CALL mma_deallocate(TVEC)
      XMAT(:,:)=XSAV(:,:)
      CALL mma_deallocate(XSAV)

      END SUBROUTINE TRACI_RPT2
