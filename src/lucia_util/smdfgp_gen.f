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
      SUBROUTINE SMDFGP_GEN(NGRP,NSMST,MXPNS,NSTFSMGP,NACTSYM,ISMDFGP)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISMDFGP(NSMST,NGRP)
      INTEGER NSTFSMGP(MXPNS,NGRP)
      INTEGER NACTSYM(NGRP)
      INTEGER IOFF
*
      NTEST = 00
*
      IF(NTEST.ge.1000) Then
        write(6,*) 'Entering  SMDFGP_GEN   '
        write(6,*) 'NGRP : ', NGRP
        write(6,*) 'NSMST: ', NSMST
        write(6,*) 'MXPNS', MXPNS
        write(6,*) 'INPUT: NSTFSMGP'
        Do ISYM = 1,NSMST
          write(6,'(40I5)') (NSTFSMGP(ISYM,IGRP),IGRP=1,NGRP)
        End Do
      End IF
*
      Do IGRP = 1,NGRP
        IOFF = 0
        NACTSYM(IGRP) = IOFF
        Do ISYM = 1,NSMST
          ISMDFGP(ISYM,IGRP) = 0
          If(NSTFSMGP(ISYM,IGRP).ne.0) Then
            IOFF = IOFF+1
            ISMDFGP(IOFF,IGRP) = ISYM
          End If
        End Do
        NACTSYM(IGRP) = IOFF
      End Do
*
      IF(NTEST.ge.100) Then
        write(6,*) 'Number of Active Symm per GRP:'
        write(6,*) (NACTSYM(IGRP),IGRP=1,NGRP)
        write(6,*) 'Symmetries allowed by each group:'
        Do ISYM = 1,NSMST
          write(6,'(40I2)') (ISMDFGP(ISYM,IGRP),IGRP=1,NGRP)
        End Do
      End IF
*
      RETURN
      END
