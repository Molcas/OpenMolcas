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
      SUBROUTINE GETDPREF(DREF,NDREF,PREF,NPREF)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, Half
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT
      use pt2_guga, only: NG1,NG2
      IMPLICIT None
      integer(kind=iwp), intent(in):: NDREF, NPREF
      real(kind=wp), intent(out):: DREF(NDREF), PREF(NPREF)

      real(kind=wp), ALLOCATABLE:: G1(:), G2(:)
      integer(kind=iwp) I,J,IJ,IJT,IJKLT,N2,JI,KLT,K,L,KL,LK,IJKL,IJLK,
     &                  JILK,JIKL
      real(kind=wp) P1, P2

* Get active 1-density and 2-density matrices GAMMA1 and
* GAMMA2, and construct DREF and PREF which are in a tringular
* storage.

C Remember: NDREF=1 if NASHT=0. Similar NPREF.
      DREF(1)=Zero
      PREF(1)=Zero
      IF(NASHT.EQ.0) RETURN
c Active density 1-matrix:
      CALL mma_allocate(G1,NG1,LABEL='G1')
      CALL PT2_GET(NG1,'GAMMA1',G1)
      DO I=1,NASHT
       DO J=1,I
        IJ=(I*(I-1))/2+J
        DREF(IJ)=G1(I+NASHT*(J-1))
       END DO
      END DO
      CALL mma_deallocate(G1)

C CONSTRUCT PREF, 2-ELECTRON DENSITY MATRIX:
      CALL mma_allocate(G2,NG2,LABEL='G2')
      CALL PT2_GET(NG2,'GAMMA2',G2)
      IJT=0
      IJKLT=0
      N2=NASHT**2
      DO I=1,NASHT
        DO J=1,I
          IJT=IJT+1
          IJ=I+NASHT*(J-1)
          JI=J+NASHT*(I-1)
          KLT=0
          DO K=1,NASHT
            DO L=1,K
              KLT=KLT+1
              IF(KLT.GT.IJT) CYCLE
              IJKLT=IJKLT+1
              KL=K+NASHT*(L-1)
              LK=L+NASHT*(K-1)

              P1=Half*G2(IJ+N2*(KL-1))
              P2=Half*G2(IJ+N2*(LK-1))
              IF(IJ.GE.KL) THEN
                IJKL=(IJ*(IJ-1))/2+KL
              ELSE
                IJKL=(KL*(KL-1))/2+IJ
              END IF
              IF(IJ.GE.LK) THEN
                IJLK=(IJ*(IJ-1))/2+LK
              ELSE
                IJLK=(LK*(LK-1))/2+IJ
              END IF
              JIKL=(JI*(JI-1))/2+KL
              JILK=(JI*(JI-1))/2+LK
              PREF(IJKL)=P1
              PREF(IJLK)=P2
              PREF(JIKL)=P2
              PREF(JILK)=P1
            END DO
          END DO
        END DO
      END DO

      CALL mma_deallocate(G2)

      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(u6,*)' GETDPREF has constructed DREF and PREF.'
       CALL XFLUSH(6)
      END IF

      END SUBROUTINE GETDPREF
