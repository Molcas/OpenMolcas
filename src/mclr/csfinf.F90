!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      Subroutine CsfInf(lSym,iSpin,MS,iSPC,iPrnt,nsym)
      use Str_Info, only: STR,CNSM,CFTP,DFTP,DTOC,NELEC,NOCTYP
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: i1,iAnders,lConf,llDET
      use MCLR_Data, only: iRefSM,iDC,MS2,PSSIGN
      use MCLR_Data, only: LuCSF2SD
      use MCLR_Data, only: IASTFI,IBSTFI,ISMOST,MNR1IC,MXR3IC,NELCI
      use MCLR_Data, only: NACOB,NORB1,NORB2,NORB3
      use MCLR_Data, only: MAXOP,MINOP,NCNATS
      use CandS, only: ICSM,ISSM,ICSPC,ISSPC
      use input_mclr, only: nIrrep
!
      Implicit None
      Integer lSym,iSpin,MS,iSPC,iPrnt,nsym
!
      integer idum(1)
      Integer, Allocatable:: SIOIO(:), SBLTP(:), IOOS1(:),              &
     &                       NOOS1(:)
      Integer NEL,IATP,IBTP,NOCTPA,NOCTPB,MNELR1,MXELR3,NOOS,IA,ISYM,   &
     &        NCOMB,LLCSF
!
!
!.... Sorry about this  but this is just to tell the program
!     that no CSF<->SD coefficents is in core
      i1=-9
      iAnders=-9
      ICSM   = lSym
      ISSM   = lSym
      ICSPC  = 1
      ISSPC  = 1
      NEL    = NELCI(ISPC)
      IATP   = IASTFI(ISPC)
      IBTP   = IBSTFI(ISPC)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      MNELR1 = MNR1IC(ISPC)
!     MXELR3 = MNR1IC(ISPC)
      MXELR3 = MXR3IC(ISPC)
      iRefSm=lsym
!.... Obtain OOS pointer array
      CALL mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,       &
     &            Str(IATP)%EL1,Str(IATP)%EL3,                          &
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,                          &
     &            SIOIO,IPRNT)
      CALL mma_allocate(SBLTP,nIrrep,Label='SBLTP')
      NOOS = NOCTPA*NOCTPB*nIrrep
      CALL mma_allocate(IOOS1,NOOS,Label='IOOS1')
      CALL mma_allocate(NOOS1,NOOS,Label='NOOS1')
      CALL INTCSF(NACOB,NEL,iSpin,MS2,                                  &
     &            NORB1,NORB2,NORB3,MNELR1,MXELR3,                      &
     &            LLCSF,1,0,PSSIGN,IPRNT,lconf,lldet)
!
!     Calculate CG COEFFICENTS ETC
!
      CALL CSDTMT(DFTP,CFTP,DTOC,PSSIGN,IPRNT)

!
!     Calculate the reordering vector and write it to disk
!

      iA=0
      Do iSym=1,nSym
!.OOS arrayy
          CALL ZBLTP(ISMOST(1,ISYM),nIrrep,IDC,SBLTP,idum)
          CALL ZOOS(ISMOST(1,ISYM),SBLTP,                               &
     &          nIrrep,SIOIO,                                           &
     &          Str(IATP)%NSTSO,Str(IBTP)%NSTSO,                        &
     &          NOCTPA,NOCTPB,idc,IOOS1,NOOS1,NCOMB,0)
!EAW           CALL CNFORD(CNSM(1)%ICTS,CNSM(1)%ICONF,
!    &                     iSym,NACOB,DFTP,
!    &          NCNATS(1,ISYM),NEL,0,0,IDUM,IDUM,
!    &          IASTFI(ISPC),IBSTFI(ISPC),IOOS1,
!    &          NORB1,NORB2,NORB3,MNELR1,MXELR3,
!    &          NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),
!    &          MINOP,MAXOP,IPRNT)
           CALL CNFORD(CNSM(1)%ICTS,CNSM(1)%ICONF,                      &
     &                 iSym,NACOB,DFTP,                                 &
     &          NCNATS(1,ISYM),NEL,0,0,IDUM,IDUM,                       &
     &          IASTFI(ISPC),IBSTFI(ISPC),IOOS1,                        &
     &          NORB1,NORB2,NORB3,MNELR1,MXELR3,                        &
     &          NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),                &
     &          MINOP,MAXOP,PSSIGN, IPRNT)
!
           CALL iDAFILE(LUCSF2SD,1,CNSM(1)%ICTS,lldet,iA)
           CALL iDAFILE(LUCSF2SD,1,CNSM(1)%ICONF,lconf,iA)
      End Do
!
!
      CALL mma_deallocate(CNSM(1)%ICTS)
      CALL mma_deallocate(CNSM(1)%ICONF)
      CALL mma_deallocate(IOOS1)
      CALL mma_deallocate(NOOS1)
      Call mma_deallocate(SBLTP)
      Call mma_deallocate(SIOIO)
!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(MS)
      END Subroutine CsfInf
