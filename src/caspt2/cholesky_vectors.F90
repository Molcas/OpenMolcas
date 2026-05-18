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
      Subroutine Cholesky_Vectors(MODE,ITK,ITQ,JSYM,Array,mArray,nArray,&
     &                            IBSTA,IBEND)

      USE CHOVEC_IO, only: NVLOC_CHOBATCH, IDLOC_CHOGROUP, NPQ_CHOTYPE
      use caspt2_module, only: NSYM
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: MODE, ITK, ITQ, JSYM, mArray,    &
     &  IBSTA, IBEND
      integer(kind=iwp), intent(out) :: nArray
      real(kind=wp), intent(_OUT_) :: Array(mArray)

      integer(kind=iwp) :: ICASE, LKETSM, LUCDER, ISYK, NQK, IB, NV,    &
     &  NKETSM, IDISK

      ! ugly hack to convert separate k/q orbital types into a specific
      ! case
      ICASE=ITK*ITQ
      IF (ICASE == 3) THEN
        ICASE=4
      ELSE
        ICASE=ICASE/2
      END IF

      LKETSM=1
      LUCDER = 63 ! tentative
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK == 0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          !! MODE = 1: write
          !! MODE = 2: read
          CALL DDAFILE(LUCDER,MODE,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
!
      End Subroutine Cholesky_Vectors
