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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE FBLOCK(FIFA,NO,NI,NA,NS,FIT,FTI,FIA,FAI,FTA,FAT)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
      INTEGER NO,NI,NA,NS
      REAL*8 FIFA((NO*(NO+1))/2)
      REAL*8 FIT(NI,NA),FIA(NI,NS),FTA(NA,NS)
      REAL*8 FTI(NA,NI),FAI(NS,NI),FAT(NS,NA)

      INTEGER IT,II,IA,ITTOT,IATOT,ITI,IAI,IAT

C Extract three rectangular submatrices FIT, FIA and FTA from the
C triangular matrix FIFA.
C SVC: add transposed submatrices to avoid complicated strides in the
C low-level sgm subroutines

      DO IT=1,NA
        ITTOT=NI+IT
        DO II=1,NI
          ITI=(ITTOT*(ITTOT-1))/2+II
          FIT(II,IT)=FIFA(ITI)
          FTI(IT,II)=FIFA(ITI)
        END DO
      END DO
      DO IA=1,NS
        IATOT=NI+NA+IA
        DO II=1,NI
          IAI=(IATOT*(IATOT-1))/2+II
          FIA(II,IA)=FIFA(IAI)
          FAI(IA,II)=FIFA(IAI)
        END DO
      END DO
      DO IA=1,NS
        IATOT=NI+NA+IA
        DO IT=1,NA
          ITTOT=NI+IT
          IAT=(IATOT*(IATOT-1))/2+ITTOT
          FTA(IT,IA)=FIFA(IAT)
          FAT(IA,IT)=FIFA(IAT)
        END DO
      END DO

      RETURN
      END
