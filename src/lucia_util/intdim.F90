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
      SUBROUTINE INTDIM(IPRNT)
!
! Number of integrals and storage mode
!
      use lucia_data, only: NINT1,I12S,I34S,I1234S,NINT2,NINT2_NO_CCSYM,&
     &                      NBINT1,NBINT2
      use lucia_data, only: PNTGRP,NSMOB
      use lucia_data, only: NTOOBS
      use lucia_data, only: MXPOBS
      use csm_data, only: ITSDX,ITSSX,NSMSX
      use csm_data, only: ADSXA,SXDXSX

      IMPLICIT NONE
      INTEGER IPRNT

      INTEGER I12,I34,I1234,NDXFSM,NSXFSM

!.1 : Number of one-electron integrals
      NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
!.2 : Number of two-electron integrals
        IF(PNTGRP.EQ.1) THEN
!. Full eightfold symmetry can be used
          I12S = 1
          I34S = 1
          I1234S = 1
!        ELSE
!. Only symmetry between 12 and 34
!          I12S = 0
!          I34S = 0
!          I1234S = 1
        END IF
        NINT2 = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,         &
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12S,I34S,I1234S,     &
     &                  IPRNT )
!. Number of integrals without complex conjugation symmetry
!. ( used for T1 transformed Hamiltonian)
       I12 = 0
       I34 = 0
       I1234 = 1
       NINT2_NO_CCSYM = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS, &
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12,I34,I1234,        &
     &                  IPRNT )
!       IF(ISIMTRH.EQ.1) THEN
!         WRITE(6,*)
!     &   ' Number of two-electron integrals in exp(-T1)Hexp(T1) ',
!     &     NINT2_NO_CCSYM
!       END IF
!. Number of symmetry blocks of one- and two-electron integrals
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
      END SUBROUTINE INTDIM
