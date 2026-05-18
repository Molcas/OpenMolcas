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
      Subroutine Get_Cholesky_Vectors(ITK,ITQ,JSYM,                     &
     &                                Array,mArray, nArray,             &
     &                                IBSTA,IBEND)
      use definitions, only: iwp, wp
      USE CHOVEC_IO, only: NPQ_CHOTYPE, NVLOC_CHOBATCH, IDLOC_CHOGROUP
      use caspt2_global, only: LUDRA
      use caspt2_module, only: NSYM
      IMPLICIT None
      integer(kind=iwp), Intent(in):: ITK,ITQ,JSYM,IBSTA,IBEND
      integer(kind=iwp), Intent(in):: mArray
      integer(kind=iwp), Intent(Out):: nArray
      real(kind=wp), intent(Out):: Array(mArray)

      integer(kind=iwp) ICASE, LKETSM, ISYK, NQK, IB, NV, NKETSM, IDISK

      ! ugly hack to convert separate k/q orbital types into a specific
      ! case
      ICASE=ITK*ITQ
      IF (ICASE.EQ.3) THEN
        ICASE=4
      ELSE
        ICASE=ICASE/2
      END IF

      LKETSM=1
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK.EQ.0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          CALL DDAFILE(LUDRA,2,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
!
      End Subroutine Get_Cholesky_Vectors
