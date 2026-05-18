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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      SUBROUTINE CnstCLag(IFF,nLev,NG3,NCONF,CLag,DG1,DG2,DG3,DF1,DF2,  &
     &                    DF3,DEPSA,G1,G2,G3)

      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: iPrGlb
      use PrintLevel, only: VERBOSE
      use sguga, only: L2ACT
      use caspt2_global, only: LUCIEX, IDTCEX, LUSOLV
      use definitions, only: wp, iwp, byte, u6
      use caspt2_module, only: STSYM, NSTATE, MSTATE, JSTATE,           &
     &                         ISCF, EPSA
      use Constants, only: One
      use caspt2_module, only: ETA, CITHR

      implicit none

      integer(kind=iwp), intent(in) :: IFF, nLev, NG3, NCONF
      real(kind=wp), intent(in) :: G1(nLev**2), G2(nLev**4), G3(NG3)
      real(kind=wp), intent(inout) :: CLag(nConf), DG1(nLev**2),        &
     &  DG2(nLev**4), DG3(NG3), DF1(nLev**2), DF2(nLev**4), DF3(NG3),   &
     &  DEPSA(nLev**2)

      integer(kind=iwp) :: ILEV, ILUID, IDCI
      integer(kind=byte), allocatable :: idxG3(:,:)
      real(kind=wp), allocatable :: CI1(:)

      real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, CPTF10, TIOE, TIOTF0,   &
     &                 TIOTF10

      IF (IFF == 1) THEN
! ORBITAL ENERGIES IN CI-COUPLING ORDER:
        DO ILEV=1,NLEV
          ETA(ILEV)=EPSA(L2ACT(ILEV))
        END DO
      END IF

!-SVC20100831: allocate local G3 matrices
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

      call mma_allocate(CI1,NCONF,LABEL='CI')
      If (ISCF == 0) Then
        if (iff == 1) then
          IDCI=IDTCEX(JSTATE)
          CALL DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
        else
!         Call LoadCI_XMS('C',1,nConf,nState,CI1,JSTATE,U0)
        end if
        IF (IPRGLB >= VERBOSE) THEN
          WRITE(u6,*)
          IF (NSTATE > 1) THEN
            WRITE(u6,'(A,I4)')                                          &
     &      ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
          ELSE
            WRITE(u6,*)' With new orbitals, the CI array is:'
          END IF
          CALL PRWF_CP2(STSYM,NCONF,CI1,CITHR)
        END IF
      Else
        CI1(1) = One
      End If

      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      If (ISCF == 0) Then
        CALL DERFG3(CI1,NCONF,NLEV,NG3,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,    &
     &              DEPSA,G1,G2)
      Else
        CALL DERSPE(NLEV,NG3,DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
      End If
      CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
      IF (IPRGLB >= VERBOSE) THEN
        CPUT =CPTF10-CPTF0
        WALLT=TIOTF10-TIOTF0
        write(u6,*)
        write(u6,'(a,2f10.2)')' DERFG3  : CPU/WALL TIME=', cput,wallt
      END IF

      call mma_deallocate(CI1)
      call mma_deallocate(idxG3)

      Return

      End Subroutine CnstCLag
