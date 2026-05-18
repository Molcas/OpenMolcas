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

      SUBROUTINE POLY1_CLagT(CI1,CI2,CLag1,CLag2,RDMEIG,Scal)
      use sguga, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: nConf
      use caspt2_module, only: MxCI, iAdr10, cLab10

      IMPLICIT NONE
! PER-AAKE MALMQUIST, 92-12-07
! THIS PROGRAM CALCULATES THE 1-EL DENSITY
! MATRIX FOR A CASSCF WAVE FUNCTION.

      real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF), RDMEIG(*),   &
     &                             Scal
      real(kind=wp), intent(inout) :: CLag1(*), CLag2(*)

      real(kind=wp),allocatable :: SGM1(:)
      integer(kind=iwp) :: nLev, I

      nLev = SGS%nLev

      IF(NLEV > 0) THEN
        call mma_allocate(SGM1,MXCI,Label='SGM1')
        CALL DENS1T_RPT2_CLag(CI1,CI2,SGM1,                             &
     &                        CLag1,CLag2,RDMEIG,Scal,nLev)
      END IF

! REINITIALIZE USE OF DMAT.
! The fields IADR10 and CLAB10 are kept in caspt2_module.F90
! CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
! HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
! ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV > 0) THEN
        call mma_deallocate(SGM1)
      END IF

      RETURN
      end subroutine POLY1_CLagT
