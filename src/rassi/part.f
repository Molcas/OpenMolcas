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
* Copyright (C) 1984,1989, Per Ake Malmqvist                           *
************************************************************************
      SUBROUTINE PART(SXY,TRA1,TRA2)
      use rasdef, only: NRS1, NRS2, NRS3
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NSXY,NTRA,NOSH,NISH

      IMPLICIT NONE
      Real*8 SXY(NSXY),TRA1(NTRA),TRA2(NTRA)
      Integer NSIZE(4)
C  PURPOSE: SXY CONTAINS THE NONSECONDARY PART OF THE MO OVERLAP
C  MATRIX. UPON RETURN, TRA1 AND TRA2 WILL CONTAIN THE COEFFICIENTS
C  FOR SEQUENTIAL SINGLE-ORBITAL TRANSFORMATIONS (VI.2, MY IJQC ARTICLE)
C  TO BIORTHONORMAL ORBITALS. SXY, TRA1 AND TRA2 ARE SYMMETRY-BLOCKED.
C  ORIGINAL VERSION, MALMQUIST 84-04-04
C  RASSCF VERSION,   MALMQUIST 89-11-15
      Real*8, allocatable:: ScrMat(:), ScrBuf(:)
      Integer, allocatable:: ScrPiv(:)
      INTEGER II,ISY,N,NBLOCK,NDIMEN,NOMAX

      NOMAX=0
      DO ISY=1,NSYM
        NOMAX=MAX(NOSH(ISY),NOMAX)
      ENDDO
      CALL mma_allocate(SCRMAT,NOMAX*NOMAX,Label='ScrMat')
      CALL mma_allocate(SCRPIV,2*NOMAX,Label='ScrPiv')
      CALL mma_allocate(SCRBUF,NOMAX,Label='ScrBuf')
      II=1
      DO ISY=1,NSYM
        NDIMEN=NOSH(ISY)
        IF(NDIMEN.EQ.0) cycle
        NBLOCK=0
        N=NISH(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS1(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS2(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS3(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        CALL PART1(NDIMEN,NBLOCK,NSIZE,SXY(II),TRA1(II),TRA2(II),
     &             SCRMAT,SCRPIV,SCRBUF)
        II=II+NDIMEN**2
      ENDDO
      CALL mma_deallocate(SCRMAT)
      CALL mma_deallocate(SCRPIV)
      CALL mma_deallocate(SCRBUF)

      END SUBROUTINE PART
