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
      SUBROUTINE PNT4DM(   NSMOB,   NSMSX,  MXPOBS,   NO1PS,   NO2PS,
     &                     NO3PS,   NO4PS,   IDXSM,   ADSXA,  SXDXSX,
     &                      IS12,    IS34,  IS1234,   IPNTR,   ISM4A,
     &                     ADASX)
*
* Pointer for 4 dimensionl array with total symmetry IDXSM
* Pointer is given as 3 dimensional array corresponding
* to the first 3 indeces
* Symmetry of last index is give by ISM4
*
* IS12 (0,1,-1)   : Permutational symmetry between indeces 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between indeces 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between indeces 12 and 34
*
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER ADASX(MXPOBS,MXPOBS)
*. Specific input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*.Output
      INTEGER IPNTR(NSMOB,NSMOB,NSMOB),ISM4A(NSMOB,NSMOB,NSMOB)
*
      CALL ISETVC(IPNTR,0,NSMOB ** 3 )
      CALL ISETVC(ISM4A,0,NSMOB ** 3 )
*
C?    WRITE(6,*) 'NO1PS NO2PS NO3PS NO4PS '
C?    CALL IWRTMA(NO1PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO2PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO3PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO4PS,1,NSMOB,1,NSMOB)
      IOFF= 1
      N12 = 0
      N34 = 0
*
      DO 10 I1SM = 1, NSMOB
        DO 20 I2SM = 1, NSMOB
          I12SM = ADASX(I1SM,I2SM)
          I34SM = SXDXSX(I12SM,IDXSM)
          IF(I34SM.EQ.0) GOTO 20
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 20
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMOB+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
          DO 30 I3SM = 1, NSMOB
            I4SM = ADSXA(I3SM,I34SM)
            IF(I4SM.EQ.0) GOTO 30
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 30
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMOB+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 30
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF+ N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF + N12*(N12+1)/2
            ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF=  IOFF+ N12*(N12-1)/2
            END IF
C?          WRITE(6,*) ' I1SM I2SM I3SM I4SM    IOFF'
C?          WRITE(6,'(1H ,4I4,I9)')   I1SM,I2SM,I3SM,I4SM,IOFF
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*
*
C?    WRITE(6,*) ' PNT4DM , 64 elemets of IPNTR '
C?    call IWRTMA(IPNTR,1,64,1,64)
      NTEST = 0
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Length of 4 index array ', IOFF - 1
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NSMSX)
      END
