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
      SUBROUTINE POLY0
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

      INTEGER I,IT,ITABS,ILEV,ISYM
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
      INTEGER IQ
#endif

      NLEV=NASHT
C ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
C PAM060612: With true RAS space, the orbitals must be ordered
C first by RAS type, then by symmetry.
      ITABS=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
! Quan: Bug in LEVEL(ITABS) and L2ACT
          if (DoCumulant) then
             do iq=1,NLEV
               LEVEL(iq)=iq
               L2ACT(iq)=iq
             enddo
          endif
#endif
          ILEV=LEVEL(ITABS)
          ISM(ILEV)=ISYM
        END DO
      END DO

C INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
      DO I=1,8
        NCSF(I)=0
      END DO
      NCSF(STSYM)=1

      IF ((.NOT.DoCumulant).AND.
     &    (NACTEL.GT.0).AND.(ISCF.EQ.0)) CALL GINIT_CP2
      MXCI=1
      DO I=1,NSYM
        MXCI=MAX(MXCI,NCSF(I))
      END DO
C NOTE: AT THIS POINT, WE HAVE ALLOCATED MEMORY SPACE FOR SGUGA USE:
C MVL,MVR,NOW,IOW,NOCP,IOCP,NOCSF,IOCSF,ICASE,ICOUP,VTAB,TMP.


      RETURN
      END
