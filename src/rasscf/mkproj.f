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
      implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "warnings.h"
#include "WrkSpc.fh"
      DIMENSION CRVEC(NCRVEC), CMO(NTOT2)
      dimension TUVX(*)

      NFI=NFRO(1)+NISH(1)
      NA=NASH(1)
      NB=NBAS(1)

      CALL GETMEM('CS_TMP','ALLO','REAL',LCS_TMP,NB)
*      write(6,*) 'MKPROJ test: Overlaps active/core :'
      DO IT=1,NA
        WORK(LCS_TMP-1+IT)=DDOT_(NB,CMO((IT-1)*NB+1),1,CRVEC,1)
*        write(6,'(1x,i5,f16.8)') it,WORK(LCS_TMP-1+IT)
      END DO
*      write(6,*)' Before shift TUVX(1):',TUVX(1)

      IPOS=0
      ITU=0
      DO IT=1,NA
        DO IU=1,IT
          CS_TU=WORK(LCS_TMP-1+IT)*WORK(LCS_TMP-1+IU)
          ITU=ITU+1
          IVX=0
          DO IV=1,IT
            CS_TUV=CS_TU*WORK(LCS_TMP-1+IV)
            IXMX=IV
            IF(IT.eq.IV) IXMX=IU
            DO IX=1,IXMX
              IVX=IVX+1
              IPOS=IPOS+1
              COREPROJ=CS_TUV*WORK(LCS_TMP-1+IX)
              TUVX(IPOS)=TUVX(IPOS)+CORESHIFT*COREPROJ
            END DO
          END DO
        END DO
      END DO

      CALL GETMEM('CS_TMP','FREE','REAL',LCS_TMP,NB)
*      write(6,*)' IPOS, NCRVEC:',IPOS, NTOT, NCRVEC
*      write(6,*)'  After shift TUVX(1):',TUVX(1)
*      call xflush(6)

      RETURN
      END
