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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

      SUBROUTINE RHSOD_H_NOSYM(IVEC)
      use definitions, only: iwp, wp
      use constants, only: Zero, Half, One, Three
      USE SUPERINDEX, only: MIGEJ, MAGEB, MIGTJ, MAGTB
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NSSHT, NAGEB, NIGEJ, NAGTB, NIGTJ

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOSYM(8,8)
      real(kind=wp), ALLOCATABLE:: CHOBUF(:)
      real(kind=wp) AJCL, ALCJ, HMACJL, HPACJL, SCL
      integer(kind=iwp) IA, IASTA, IAEND, IISTA, IIEND, MW, IAGEB,      &
     &                  IAGTB, IC, iCASE, IDX, IJ, IJGEL, IJGTL, IL,    &
     &                  lg_W, LIJOFF, LILOFF, NAS, NBLOCK, NCHOBUF,     &
     &                  NIS, NV, NW
!      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), PARAMETER :: NOSYM = 1
      real(kind=wp), ALLOCATABLE :: AIBJ(:,:)
      real(kind=wp), parameter:: SQRT3=SQRT(Three), SQRTH=SQRT(Half)

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case H'
      END IF

!***********************************************************************
! Case H:
!   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
!   WM(jl,ac)=((ajcl)-(alcj))*SQRT(Three)
!***********************************************************************

      NV=NVTOT_CHOSYM(NOSYM)
      Call mma_ALLOCATE(AIBJ,NSSHT,NSSHT,Label='AIBJ')
      NBLOCK=NV*NSSHT

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
      CALL CHOVEC_SIZE(4,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,Label='CHOBUF')

      CALL CHOVEC_READ(4,CHOBUF,NCHOBUF)

      iCASE=12

      NAS=NAGEB(NOSYM)
      NIS=NIGEJ(NOSYM)
      NW=NAS*NIS

      IF(NW/=0) Then

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGEL=IISTA,IIEND
        IJ=MIGEJ(1,IJGEL)
        IL=MIGEJ(2,IJGEL)
        LIJOFF=1+NBLOCK*(IJ-1)
        LILOFF=1+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,                             &
     &              One,CHOBUF(LIJOFF),NV,CHOBUF(LILOFF),NV,            &
     &              Zero,AIBJ,NSSHT)
        DO IAGEB=IASTA,IAEND ! these are always all elements
          IA=MAGEB(1,IAGEB)
          IC=MAGEB(2,IAGEB)
          AJCL=AIBJ(IA,IC)
          ALCJ=AIBJ(IC,IA)
! HP(ac,jl)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
          SCL=One
          IF (IA==IC) SCL=SCL*SQRTH
          IF (IL==IJ) SCL=SCL*SQRTH
          HPACJL=SCL*(AJCL+ALCJ)
! write element HP(ac,jl)
          IDX=IAGEB+NAS*(IJGEL-IISTA)
#ifdef _MOLCAS_MPP_
          DBL_MB(MW+IDX-1)=HPACJL
#else
          GA_Arrays(lg_w)%A(IDX)=HPACJL
#endif
        END DO
      END DO
!***********************************************************************
      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (lg_W)
!***********************************************************************
      End If

      iCASE=13

      NAS=NAGTB(NOSYM)
      NIS=NIGTJ(NOSYM)
      NW=NAS*NIS

      IF(NW/=0) Then

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGTL=IISTA,IIEND
        IJ=MIGTJ(1,IJGTL)
        IL=MIGTJ(2,IJGTL)
        LIJOFF=1+NBLOCK*(IJ-1)
        LILOFF=1+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,                             &
     &              One,CHOBUF(LIJOFF),NV,CHOBUF(LILOFF),NV,            &
     &              Zero,AIBJ,NSSHT)
        DO IAGTB=IASTA,IAEND ! these are always all elements
          IA=MAGTB(1,IAGTB)
          IC=MAGTB(2,IAGTB)
          AJCL=AIBJ(IA,IC)
          ALCJ=AIBJ(IC,IA)
! HM(ac,jl)=((ajcl)-(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
          SCL=SQRT3
          HMACJL=SCL*(AJCL-ALCJ)
! write element HM(ac,jl)
          IDX=IAGTB+NAS*(IJGTL-IISTA)
#ifdef _MOLCAS_MPP_
          DBL_MB(MW+IDX-1)=HMACJL
#else
          GA_Arrays(lg_W)%A(IDX)=HMACJL
#endif
        END DO
      END DO
!***********************************************************************

      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (lg_W)
      End If
!***********************************************************************

      CALL mma_deallocate(CHOBUF)

      call mma_DEALLOCATE(AIBJ)

      END SUBROUTINE RHSOD_H_NOSYM
