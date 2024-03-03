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
      SUBROUTINE MKISM(SGS)
      use UnixInfo, only: ProgName
      use gugx, only: LEVEL
      use stdalloc, only: mma_allocate
      use Struct, only: SGStruct
      Implicit None
      Type(SGStruct) SGS
#include "rassi.fh"
      Integer ITABS, ISYM, IT, ILEV, nSym

      If (ProgName(1:6)=='rassi') Then

         nSym=SGS%nSym
         SGS%NLEV=NASHT ! Total number of active orbitals
! Allocate Level to Symmetry table ISm:
         Call mma_allocate(SGS%ISm,SGS%nLev,Label='SGS%ISm')
         ITABS=0
         DO ISYM=1,NSYM
           DO IT=1,NASH(ISYM)
             ITABS=ITABS+1
             ILEV=LEVEL(ITABS)
             SGS%ISM(ILEV)=ISYM
           END DO
         END DO
      Else

         Call mkNSM()

      End If

      END SUBROUTINE MKISM

      SUBROUTINE MKNSM()
!     PUPROSE: CREATE THE SYMMETRY INDEX VECTOR
!
      use gugx, only: SGS
      use stdalloc, only: mma_allocate
      IMPLICIT None
!
! to get some dimensions
#include "rasdim.fh"
! NSM form rasscf,fh
#include "rasscf.fh"
! NSYM from general.fh
#include "general.fh"
! NGAS and NGSSH from gas.fh
#include "gas.fh"
!
      Integer IGAS, ISYM, LEV, NLEV, NSTA

      NLEV=0
      DO IGAS=1,NGAS
        DO ISYM=1,NSYM
          NSTA=NLEV+1
          NLEV=NLEV+NGSSH(IGAS,ISYM)
          DO LEV=NSTA,NLEV
            NSM(LEV)=ISYM
          END DO
        END DO
      END DO

      If (SGS%nSym/=0) Then
         SGS%nLev=nLev
         Call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')
         SGS%ISM(1:nLev)=NSM(1:nLev)
      End If

      END SUBROUTINE MKNSM
