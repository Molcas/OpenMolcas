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
* Copyright (C) 2012, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE NEXT_SYM_DISTR_NEW(NSMST,NGRP,KGRP,NGAS,
     &                              ISYM,ISYM_TOT,IFIRST,
     &                              NONEW,ISMDFGP,NACTSYM,ISMSCR)
*
* Giovanni Li Manni. Geneva, February 2012
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER ISMDFGP(NSMST,NGRP),NACTSYM(NGRP),ISMSCR(NGRP)
      INTEGER KGRP(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*
      NTEST = 00
*
      IF(NTEST.ge.100) then
        write(6,*) '***************************'
        write(6,*) 'INPUT IN NEXT_SYM_DISTR_NEW'
        write(6,*) '***************************'
        write(6,*) 'NGRP :', NGRP
        write(6,*) 'KGRP :'
        write(6,'(40I3)') (KGRP(i),i=1,NGAS)
        write(6,*) 'NACTSYM(NGRP) :'
        write(6,'(40I2)') (NACTSYM(i),i = 1,NGRP)
        write(6,*) 'ISMDFGP(ISYM,NGRP) :'
        do Ism = 1, NSMST
          write(6,'(40I2)') (ISMDFGP(ism,i),i = 1, NGRP)
        end do
        write(6,*) 'IFIRST', IFIRST
      End IF
*
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS
          ISMSCR(IGAS) = 1
          ISYM(IGAS) = ISMDFGP(1,KGRP(IGAS))
        END DO
        NONEW = 0
      END IF
*
      IF(NTEST.ge.100) then
        write(6,*) 'Symmetry distribution :'
        write(6,'(40I2)') (ISYM(IGAS),IGAS=1,NGAS)
      End IF
*
1001  CONTINUE
      IF(IFIRST.EQ.0) then
        CALL NXTDIST(NSMST,NGRP,NGAS,KGRP,ISMDFGP,ISMSCR,NACTSYM,NONEW)
        Do IGAS = 1, NGAS
          ISYM(IGAS) = ISMDFGP(ISMSCR(IGAS),KGRP(IGAS))
        End Do
      END if
      IFIRST = 0
      IF(NONEW.EQ.0) THEN
        JSYM = ISYMSTR(ISYM,NGAS)
        IF(JSYM.ne.ISYM_TOT)GOTO 1001
      END IF
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
