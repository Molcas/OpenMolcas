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
      SUBROUTINE TRACID(       T,   LUCIN,  LUCOUT,   LUSC1,   LUSC2,
     &                     LUSC3,    VEC1,    VEC2)
*
* Transform CI vector on LUCIN with T matrix after
* Docent Malmquist's recipe. Place result as next vector on LUOUT
*
* The transformation is done as a sequence of one-electron transformations
*
* with each orbital transformation being
*
* Sum(k=0,2) ( 1/k! sum(n'.ne.n) S(n'n) E_{n'n} ) Tnn^N_n
*
* with Sn'n = T(n'n)/Tnn
*
* each transformation is
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "glbbas.fh"
#include "oper.fh"
#include "intform.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "io_util.fh"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPRDD
*. Input
      DIMENSION T(*)
*. Scratch blocks ( two of them)
      DIMENSION VEC1(*),VEC2(*)
*
      IDUM=0
*
      NTEST = 000
      LBLK = -1
      IDUM = 1
*. Transfer vector on LUCIN to LUSC1
C           COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL  COPVCD(LUCIN,LUSC1,VEC1,1,LBLK)
*. A bit of info for the sigma routine
      I_RES_AB = 0
*. Do the one-electron update
        I12 = 1
*. With 1-electron integrals in complete block form
        IH1FORM = 2
*. Transform each orbital separately
      DO K = 1, NTOOB
*. Place (T(P,K)/S(K,K)   in one-electron integral list
C                       T_ROW_TO_H(T,H,K)
        CALL T_ROW_TO_H(T,WORK(KINT1),K,TKK)
*. T_{kk}^Nk
C            T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
        CALL T_TO_NK_VEC(     TKK,       K,    ISSM,   ISSPC,   LUSC1,
     &                      LUSC2,    VEC1)
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' output from T_TO_NK'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
*. For each orbital calculate (1+T+1/2 T^2)|0>
* + T
        CALL MV7(VEC1,VEC2,LUSC1,LUSC2)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector'
          CALL WRTVCD(VEC1,LUSC2,1,LBLK)
        END IF
        ONE = 1.0D0
        CALL VECSMDP(     VEC1,     VEC2,      ONE,      ONE,    LUSC1,
     &                   LUSC2,    LUSC3,        1,     LBLK)
        CALL COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated vector'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
*. + 1/2 T^2
        CALL MV7(VEC1,VEC2,LUSC2,LUSC3)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector'
          CALL WRTVCD(VEC1,LUSC3,1,LBLK)
        END IF
        ONE = 1.0D0
        HALF  = 0.5D0
        CALL VECSMDP(     VEC1,     VEC2,      ONE,     HALF,    LUSC1,
     &                   LUSC3,    LUSC2,        1,     LBLK)
*. and transfer back to LUSC1
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated vector'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
      END DO
*. And transfer to LUCOUT
      CNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
      IF (NTEST .GT. 0) WRITE(6,*) ' Norm of transformed vector', CNORM
C?    WRITE(6,*) ' Transformed vector'
C?    CALL WRTVCD(VEC1,LUSC1,1,LBLK)
      IDISK(LUSC1)=0
C?    WRITE(6,*) ' LUCOUT LUSC1 = ', LUCOUT,LUSC1
      CALL COPVCD(LUSC1,LUCOUT,VEC1,0,LBLK)
*
*
      RETURN
      END
