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
      use Struct, only: SGStruct
      Implicit None
      Type(SGStruct) SGS

      If (ProgName(1:6)=='rassi') Then

         Call mkISm_Rassi(SGS)

      Else If (ProgName(1:4)=='mclr') Then

         Call mkISm_mclr(SGS)

      Else If (ProgName(1:6)=='caspt2') Then

         Call mkISm_cp2(SGS)

      Else

         Call mkNSM()

      End If

      END SUBROUTINE MKISM

      SUBROUTINE MKISM_MCLR(SGS)
      use Struct, only: SGStruct
      use stdalloc, only: mma_allocate
      Implicit None
      Type(SGStruct) SGS

#include "Input.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Integer :: iOrb, iSym, iBas

      SGS%NLEV=ntASh

      Call mma_allocate(SGS%ISM,SGS%nLev,Label='SGS%ISM')

      iOrb=0
      Do iSym=1,nSym
         Do iBas=1,nRs1(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs2(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs3(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do

      End SUBROUTINE MKISM_MCLR

      SUBROUTINE MKISM_RASSI(SGS)
      use gugx, only: LEVEL
      use stdalloc, only: mma_allocate
      use Struct, only: SGStruct
      Implicit None
      Type(SGStruct) SGS
#include "rassi.fh"
      Integer ITABS, ISYM, IT, ILEV, nSym


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

      END SUBROUTINE MKISM_RASSI

      SUBROUTINE mkism_cp2(SGS)

      use fciqmc_interface, only: DoFCIQMC
      use stdalloc, only: mma_allocate
      use gugx, only: L2ACT, LEVEL
      use Struct, only: SGStruct

      IMPLICIT NONE
      Type(SGStruct) SGS

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"
      Integer nLev

      INTEGER IT,ITABS,ILEV,ISYM, iq

      NLEV=NASHT
      SGS%nLev = NLEV
      Call mma_allocate(SGS%ISM,NLEV,Label='ISM')
! ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
! PAM060612: With true RAS space, the orbitals must be ordered
! first by RAS type, then by symmetry.
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

      END SUBROUTINE mkism_cp2

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
