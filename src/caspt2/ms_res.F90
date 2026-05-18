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

      Subroutine MS_Res(MODE,IST,JST,Scal)
      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV, only: IVECC, IVECC2, IVECW
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: STSYM, NCONF, NASHT, ISCF, NSTATE
      use caspt2_module, only: MXCI
      use Constants, only: Zero
      use definitions, only: wp, iwp
!
!     Compute the derivative of E^PT2 with respct to the T amplitude
!
      implicit none

      integer(kind=iwp), intent(in) :: MODE, IST, JST
      real(kind=wp), intent(in) :: Scal

      integer(kind=iwp) :: I, NTG1, NTG2, NTG3, IDCI
      real(kind=wp) :: DVALUE, OVL, DUMMY(1)

      real(kind=wp),allocatable :: TG1(:),TG2(:),TG3(:),CI1(:),CI2(:)

! We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
! Note: Need proper allocation even if unused, sinced allocated
! arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      call mma_allocate(TG1,NTG1,Label='TG1')
      call mma_allocate(TG2,NTG2,Label='TG2')
      call mma_allocate(TG3,NTG3,Label='TG3')
      TG1(:) = Zero
      TG2(:) = Zero
      TG3(:) = Zero

      call mma_allocate(CI1,MXCI,Label='MCCI1')
      call mma_allocate(CI2,MXCI,Label='MCCI2')
      IF(ISCF == 0) THEN
! Read root vectors nr. IST and JST from LUCI.
        DO I=1,NSTATE
          IDCI=IDTCEX(I)
          IF(I == IST) THEN
            CALL DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
            IF(I == JST) THEN
              CI2(1:NCONF) = CI1(1:NCONF)
            END IF
          ELSE IF(I == JST) THEN
            CALL DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          END IF
        END DO
      END IF

      CALL MKTG3(STSYM,STSYM,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)

!! Do similar to RHS_STRANS. Multiply the solution vector (T) with
!! the overlap-like term constructe with transition density
!! matrices. The output is IVECC (MODE=1) or IVECC2 (MODE=2)
      IF (MODE == 1) THEN
        CALL MS_STRANS(IVECW,IVECC,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,   &
     &                 SCAL)
      ELSE IF (MODE == 2) THEN
        CALL MS_STRANS(IVECC,IVECC2,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,  &
     &                 SCAL)
      ELSE IF (MODE == 3) THEN
        CALL MS_STRANS(IVECW,IVECC,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,   &
     &                 SCAL)
      END IF

      call mma_deallocate(TG1)
      call mma_deallocate(TG2)
      call mma_deallocate(TG3)

      End Subroutine MS_Res
