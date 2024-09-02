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
* Copyright (C) 2019, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE MKPROJ(CRVEC,CMO,TUVX)
      use stdalloc, only: mma_allocate, mma_deallocate
      implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "warnings.h"
      Real*8 CRVEC(NCRVEC), CMO(NTOT2)
      Real*8 TUVX(*)
      Real*8, Allocatable:: CS_TMP(:)

      NA=NASH(1)
      NB=NBAS(1)

      CALL mma_allocate(CS_TMP,NB,Label='CS_TMP')
*      write(6,*) 'MKPROJ test: Overlaps active/core :'
      DO IT=1,NA
        CS_TMP(IT)=DDOT_(NB,CMO((IT-1)*NB+1),1,CRVEC,1)
*        write(6,'(1x,i5,f16.8)') it,CS_TMP(IT)
      END DO
*      write(6,*)' Before shift TUVX(1):',TUVX(1)

      IPOS=0
      ITU=0
      DO IT=1,NA
        DO IU=1,IT
          CS_TU=CS_TMP(IT)*CS_TMP(IU)
          ITU=ITU+1
          IVX=0
          DO IV=1,IT
            CS_TUV=CS_TU*CS_TMP(IV)
            IXMX=IV
            IF(IT.eq.IV) IXMX=IU
            DO IX=1,IXMX
              IVX=IVX+1
              IPOS=IPOS+1
              COREPROJ=CS_TUV*CS_TMP(IX)
              TUVX(IPOS)=TUVX(IPOS)+CORESHIFT*COREPROJ
            END DO
          END DO
        END DO
      END DO

      CALL mma_deallocate(CS_TMP)
*      write(6,*)' IPOS, NCRVEC:',IPOS, NTOT, NCRVEC
*      write(6,*)'  After shift TUVX(1):',TUVX(1)
*      call xflush(6)

      END SUBROUTINE MKPROJ
