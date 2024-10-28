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
      SUBROUTINE POLY0()

      use fciqmc_interface, only: DoFCIQMC
      use stdalloc, only: mma_allocate
      use gugx, only: SGS, L2ACT, LEVEL, CIS, EXS

      IMPLICIT NONE

#include "caspt2.fh"
#include "pt2_guga.fh"

      Integer nLev

      INTEGER I,IT,ITABS,ILEV,ISYM, iq

      if ((.NOT.DoCumulant) .and. (nactel.gt.0) .and. (iscf.eq.0)
     &      .and. (.not. DoFCIQMC)) Then
!     if ((.NOT.DoCumulant) .and. (nactel.gt.0)) Then

         call sginit_cp2(nSym,iSpin,nActEl,nHole1,nEle3,
     &                   nRas1T,nRas2T,nRas3T,SGS,CIS,EXS)

      else

      NLEV=NASHT
      SGS%nLev = NLEV
      Call mma_allocate(SGS%ISM,NLEV,Label='ISM')
C ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
C PAM060612: With true RAS space, the orbitals must be ordered
C first by RAS type, then by symmetry.
      ITABS=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
! Quan: Bug in LEVEL(ITABS) and L2ACT
          if (DoCumulant .or. DoFCIQMC) then
             do iq=1,NLEV
               LEVEL(iq)=iq
               L2ACT(iq)=iq
             enddo
          endif
          ILEV=LEVEL(ITABS)
          SGS%ISM(ILEV)=ISYM
        END DO
      END DO
C INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
         Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')
         CIS%NCSF(:)=0
         CIS%NCSF(STSYM)=1
         Call mma_allocate(EXS%ICoup,[1,3],[1,1],Label='EXS%ICoup')
         Call mma_allocate(EXS%VTab,[1,1],Label='EXS%VTab')
      endif

      MXCI=1
      DO I=1,NSYM
        MXCI=MAX(MXCI,CIS%NCSF(I))
      END DO

C NOTE: AT THIS POINT, WE HAVE ALLOCATED MEMORY SPACE FOR SGUGA USE:
C MVL,MVR,NOW,IOW,NOCP,IOCP,NOCSF,IOCSF,ICASE,ICOUP,VTAB,TMP.


      END SUBROUTINE POLY0
