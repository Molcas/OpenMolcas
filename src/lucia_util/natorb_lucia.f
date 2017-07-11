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
* Copyright (C) 1994,1998, Jeppe Olsen                                 *
*               2017, Roland Lindh                                     *
************************************************************************
      SUBROUTINE NATORB_LUCIA(   RHO1,  NSMOB, NTOPSM, NACPSM, NINPSM,
     &                          ISTOB,   XNAT, RHO1SM, OCCNUM,  NACOB,
     &                            SCR, IPRDEN)
*
* Obtain natural orbitals in symmetry blocks
*
* Jeppe Olsen, June 1994
*              Modification, Oct 94
*              Last modification, Feb. 1998 (reorder deg eigenvalues)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION RHO1(NACOB,NACOB)
      DIMENSION ISTOB(*)
      DIMENSION NTOPSM(NSMOB), NACPSM(NSMOB),NINPSM(NSMOB)
*. Output
      DIMENSION RHO1SM(*),OCCNUM(*),XNAT(*)
*. Scratch ( Largest symmetry block )
      DIMENSION SCR(*)
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDEN)
*. To get rid of annoying and incorrect compiler warnings
      IOBOFF = 0
      IMTOFF = 0
*. IOBOFF : Offset for active orbitals in symmetry order
      DO ISMOB = 1, NSMOB
        LOB = NACPSM(ISMOB)
        IF(ISMOB.EQ.1) THEN
          IOBOFF = NINPSM(1)+1
          IMTOFF = 1
        ELSE
          IOBOFF =
     &    IOBOFF + NTOPSM(ISMOB-1)-NINPSM(ISMOB-1)+NINPSM(ISMOB)
          IMTOFF = IMTOFF + NACPSM(ISMOB-1)**2
        END IF
        LOB = NACPSM(ISMOB)
*
*. Extract symmetry block of density matrix
*
        DO IOB = IOBOFF,IOBOFF + LOB-1
           DO JOB = IOBOFF,IOBOFF + LOB-1
*. Corresponding type indeces
             IOBP = ISTOB(IOB)
             JOBP = ISTOB(JOB)
             RHO1SM(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1)
     &     = RHO1(IOBP,JOBP)
           END DO
        END DO
*
        IF(NTEST.GE.2 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' Density matrix for symmetry  = ', ISMOB
          WRITE(6,*) ' ======================================='
          WRITE(6,*)
          CALL WRTMAT(RHO1SM(IMTOFF),LOB,LOB,LOB,LOB)
        END IF
*. Pack and diagonalize
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
        CALL TRIPAK(RHO1SM(IMTOFF),SCR,1,LOB,LOB)
        ONEM = -1.0D0
*. scale with -1 to get highest occupation numbers as first eigenvectors
        CALL SCALVE(SCR,ONEM,LOB*(LOB+1)/2)
        Call DCopy_(LOB**2,0.0D0,0,R,1)
        Call DCopy_(LOB,1.0D0,0,R,1+LOB)
        Call NIDiag(SCR,XNAT(IMTOFF),LOB,LOB,0)
        Call JACORD(SCR,XNAT(IMTOFF),LOB,LOB)
*
        DO  I = 1, LOB
          OCCNUM(IOBOFF-1+I) = - SCR(I*(I+1)/2)
        END DO
*. Order the degenerate eigenvalues so diagonal terms are maximized
        TESTY = 1.0D-11
        DO IOB = 2, LOB
          IF(ABS(OCCNUM(IOBOFF-1+IOB)-OCCNUM(IOBOFF-2+IOB))
     &       .LE.TESTY) THEN
            XII   = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB  ))
            XI1I1 = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB-1))
            XII1  = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB  ))
            XI1I  = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB-1))
*
            IF( XI1I.GT.XII.AND.XII1.GT.XI1I1 ) THEN
*. interchange orbital IOB and IOB -1
              CALL SWAPVE(XNAT(IMTOFF+(IOB-1)*LOB),
     &                    XNAT(IMTOFF+(IOB-1-1)*LOB),LOB)
              SS = OCCNUM(IOBOFF-1+IOB-1)
              OCCNUM(IOBOFF-1+IOB-1) = OCCNUM(IOBOFF-1+IOB)
              OCCNUM(IOBOFF-1+IOB)   = SS
              IF (NTEST .GE. 1)
     &           write(6,*) ' Orbitals interchanged ',
     &                      IOBOFF-1+IOB,IOBOFF-2+IOB
            END IF
          END IF
        END DO
*
        IF(NTEST.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*)
     &    ' Natural occupation numbers for symmetry = ', ISMOB
          WRITE(6,*)
     &    ' ==================================================='
          WRITE(6,*)
          CALL WRTMAT(OCCNUM(IOBOFF),1,LOB,1,LOB)
          IF(NTEST.GE.2 ) THEN
            WRITE(6,*)
            WRITE(6,*) ' Corresponding Eigenvectors '
            WRITE(6,*)
            CALL WRTMAT(XNAT(IMTOFF),LOB,LOB,LOB,LOB)
          END IF
        END IF
      END DO
*. ( End of loop over orbital symmetries )
*
      RETURN
      END
*
CADDE
