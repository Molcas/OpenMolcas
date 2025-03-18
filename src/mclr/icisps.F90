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
      SUBROUTINE ICISPS(IPRNT)
      Use Str_Info, only: STR, NOCTYP
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: IDC
      use MCLR_Data, only: IASTFI,IBSTFI,ISMOST,MNR1IC,MXR3IC,IACTI,    &
     &                       MNR3IC,MXR1IC,NAELCI,NBELCI,XISPSM,MXSB,   &
     &                       MXSOOB,NICISP
      use DetDim, only: MXPCSM
      use Constants, only: Zero
      use input_mclr, only: nIrrep
!
! Number of dets and combinations
! per symmetry for each type of internal space
!
! Jeppe Olsen, Winter 1991
! Last revision April 1991
      IMPLICIT None
      Integer IPRNT

! local variables
      Integer MXCEXP,ICI,ISYM,IATP,IBTP,IIDC,NTEST,MX,MXS,MXSOO,NCOMB
      Real*8 XNCOMB
!
      Integer, Allocatable:: LBLTP(:), LCVST(:)
! ====================
!. Output  XISPSM is calculated
! ====================
!
!
!
!.Local memory
      Call mma_allocate(LBLTP,nIrrep,Label='LBLTP')
!     IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
!    &Call mma_allocate(LCVST,nIrrep,Label='LCVST')
      Call mma_allocate(LCVST,nIrrep,Label='LCVST')

!. Obtain array giving symmetry of sigma v reflection times string
!. symmetry.
!     IF(IDC.EQ.3.OR.IDC.EQ.4) CALL SIGVST(LCVST,nIrrep)

!. Array defining symmetry combinations of internal strings
!. Number of internal dets for each symmetry
!            SMOST(nIrrep,nIrrep,MXPCSM,ISMOST)
        CALL SMOST_MCLR(nIrrep,nIrrep,MXPCSM,ISMOST)

      MXSB = 0
      MXSOOB = 0
      MXCEXP = 0
      XISPSM(:,:) = Zero
      DO 100 ICI = 1, NICISP
      DO  50 ISYM = 1, nIrrep
        IATP = IASTFI(ICI)
        IBTP = IBSTFI(ICI)
!MS        write(6,*) ' NRASDT : ICI IATP IBTP ',ICI,IATP,IBTP
        IF(NAELCI(ICI).EQ.NBELCI(ICI)) THEN
          IIDC = IDC
        ELSE
          IIDC = 1
        END IF
        IF(IACTI(ICI).EQ.1) THEN
          CALL ZBLTP(ISMOST(1,ISYM),nIrrep,IIDC,LBLTP,LCVST)
          CALL NRASDT(MNR1IC(ICI),MXR1IC(ICI),MNR3IC(ICI),MXR3IC(ICI),  &
     &         ISYM,nIrrep,NOCTYP(IATP),NOCTYP(IBTP),Str(IATP)%EL1,     &
     &         Str(IBTP)%EL1,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,           &
     &         Str(IATP)%EL3,Str(IBTP)%EL3,                             &
     &         NCOMB,XNCOMB,MXS,MXSOO,LBLTP)
          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB = MAX(MXSOOB,MXSOO)
          MXSB = MAX(MXSB,MXS)
          MXCEXP = MAX(MXCEXP,NCOMB)
!       ELSE
!         XISPSM(ISYM,ICI) = Zero
        END IF
   50 CONTINUE
  100 CONTINUE
      Call mma_deallocate(LCVST)
      Call mma_deallocate(LBLTP)
!
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      If (ntest.ne.0) Then
      WRITE(6,*)                                                        &
     &' Number of internal combinations per symmetry '
      WRITE(6,*)                                                        &
     & ' =========================================== '
      IF(NTEST.EQ.0) THEN
        MX = 1
      ELSE
        MX = NICISP
      END IF
!
      DO 200 ICI = 1, MX
        IF(IACTI(ICI).EQ.1) THEN
          WRITE(6,*) ' Internal CI space ', ICI
          CALL WRTMAT(XISPSM(1,ICI),1,nIrrep,1,nIrrep)
        END IF
  200 CONTINUE
      WRITE(6,*) ' Largest CI space                 ',MXCEXP
      WRITE(6,*) ' Largest symmetry block           ',MXSB
      WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
      End If
!
      END SUBROUTINE ICISPS
