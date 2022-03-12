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
! Here we construct the density matrix given the orbital
! coefficients.
      SUBROUTINE DENSI_MO(DENS,ORBCO,IS,IA,NBAS,IDIM)
      IMPLICIT Real*8 (A-H,O-Z)
      DIMENSION DENS(*),ORBCO(IDIM,*)
      IJ=0
      DO 8 I=1,NBAS
        DO 9 J=1,I
          IJ=IJ+1
          DENS(IJ)=0.0d0
   9    CONTINUE
   8  CONTINUE
      DO 10 I=IS,IS+IA-1
        IJ=0
        DO 11 J=1,NBAS
          DO 12 K=1,J
            IJ=IJ+1
            DENS(IJ)=DENS(IJ)+4.d0*ORBCO(J,I)*ORBCO(K,I)
 12       CONTINUE
          DENS(IJ)=DENS(IJ)-ORBCO(J,I)*ORBCO(J,I)*2.d0
 11     CONTINUE
  10  CONTINUE
      RETURN
      END
