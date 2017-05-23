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
      SUBROUTINE PREPSD(WFTP,TRORB,ISGSTR,ICISTR,IXSTR,LSYM,TRA,
     &                  ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,
     &                  NCONF,CI,DET)
      IMPLICIT NONE
      LOGICAL TRORB
      INTEGER ISGSTR(*),ICISTR(*),IXSTR(*)
      INTEGER ICNFTAB(*),ISPNTAB(*),ISSTAB(*),IFSBTAB(*)
      INTEGER LSYM,NCONF
      REAL*8 TRA(*),CI(*),DET(*)
      INTEGER IMODE,LCTMP
      CHARACTER*8 WFTP
#include "WrkSpc.fh"
C Purpose: Given a RASSCF wave function in Split-GUGA format
C and an orbital transformation matrix for the purpose of
C getting biorthonormal orbitals, prepare a wave function
C in the general SD format, using transformed orbitals.

      IF(TRORB) THEN
C Transform to bion basis, Split-Guga format
        CALL CITRA (WFTP,ISGSTR,ICISTR,IXSTR,LSYM,TRA,NCONF,CI)
      END IF

      IF(WFTP.EQ.'GENERAL ') THEN
C Transform SGUGA to SymmG:
        CALL GETMEM('PREPSD','ALLO','REAL',LCTMP,NCONF)
        IMODE=1
        CALL SYG2SGU(IMODE,ISGSTR,ICISTR,LSYM,ICNFTAB,ISPNTAB,
     &                  CI,WORK(LCTMP))
C Transform SymmG to Slater Dets:
        CALL SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,WORK(LCTMP),DET)
        CALL GETMEM('PREPSD','FREE','REAL',LCTMP,NCONF)
      ELSE
        DET(1)=CI(1)
      END IF
      RETURN
      END
