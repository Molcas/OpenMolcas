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
      Subroutine DerHEff(nConf,nRoots,nState,CLag,VECROT)

      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV, only: IVECW, IVECC
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: STSYM, NASHT, ISCF, JSTATE
      use caspt2_module, only: MXCI
      use Constants, only: Zero

      implicit none

      integer(kind=iwp), intent(in) :: nConf, nRoots, nState
      real(kind=wp), intent(inout) :: CLag(nConf,nRoots)
      real(kind=wp), intent(in) :: VECROT(nState)

      integer(kind=iwp) ::  IST, JST, I, NTG1, NTG2, NTG3, IDCI
      real(kind=wp) :: OVL, DUMMY(1)

      real(kind=wp),allocatable :: DTG1(:),DTG2(:),DTG3(:),CI1(:),      &
     &  CI2(:),CI3(:)
!     return

! We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
! Note: Need proper allocation even if unused, sinced allocated
! arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      call mma_allocate(DTG1,NTG1,Label='DTG1')
      call mma_allocate(DTG2,NTG2,Label='DTG2')
      call mma_allocate(DTG3,NTG3,Label='DTG3')
      DTG1(:) = Zero
      DTG2(:) = Zero
      DTG3(:) = Zero

      !! OVL will contain the derivative contribution?
      !! It should be ignored
      OVL = Zero
      CALL DerHeffX(IVECW,IVECC,NASHT,NTG3,OVL,DTG1,DTG2,DTG3)

      call mma_allocate(CI1,MXCI,Label='MCCI1')
      call mma_allocate(CI2,MXCI,Label='MCCI2')
      call mma_allocate(CI3,MXCI,Label='MCCI3')

      IF(ISCF == 0) THEN
        JST = jState
        CI1(1:NCONF) = Zero
        DO I=1,NSTATE
          IDCI=IDTCEX(I)
          IST = I
          IF (IST == JST) THEN
            CALL DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
          Else If (ABS(VECROT(IST)) <= 1.0e-12_wp) Then
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,2,CI3,NCONF,IDCI)
            CI1(1:NCONF) = CI1(1:NCONF) + VECROT(IST)*CI3(1:NCONF)
          END IF
        END DO

        CI3(1:NCONF) = Zero
        CALL DERTG3(.TRUE.,STSYM,STSYM,NCONF,NASHT,CI1,CI2,OVL,         &
     &              DTG1,DTG2,NTG3,DTG3,CI3,CLag(1,JST))

        DO I=1,NSTATE
          IST = I
          IF (IST == JST) THEN
            CYCLE
          Else If (ABS(VECROT(IST)) <= 1.0e-12_wp) Then
            CYCLE
          ELSE
            CLag(1:NCONF,IST) = CLag(1:NCONF,IST)                       &
     &        + VECROT(IST)*CI3(1:NCONF)
          END IF
        END DO
      END IF

      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
      call mma_deallocate(CI3)

      call mma_deallocate(DTG1)
      call mma_deallocate(DTG2)
      call mma_deallocate(DTG3)

      End Subroutine DerHEff
