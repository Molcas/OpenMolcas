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
      SUBROUTINE HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
!
! RASSCF program: version IBM-3090: SX section
!
! Purpose: Calculation of the Davidson Hamiltonian.
! The sigma vector HC is computed in SIGVEC
!
! ********** IBM-3090 Release 88 09 08 **********
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: DIA, SXN, BM, F1, F2, SXG, SXH, NLX
      use rasscf_global, only: NSXS

      IMPLICIT None
      Integer NDIM, NDIMH, NTRIAL
      REAL*8 C(*),HC(*),HH(*)
      REAL*8 HD(NDIM)
      Real*8, Allocatable:: XX(:), XC(:)
      Integer IST, IJ, L1, I, JST, J
      REAL*8, External :: DDot_
!
!
! COMPUTE SIGMA VECTORS
!
      IST=1+NDIMH*NDIM

      CALL mma_allocate(XX,NLX,Label='XX')
      CALL mma_allocate(XC,NSXS,Label='XC')

      CALL SIGVEC(C(IST),HC(IST),HD,BM,SXN,SXG,SXH,DIA,                 &
     &            F1,F2,XX,XC,NTRIAL)

      CALL mma_deallocate(XX)
      CALL mma_deallocate(XC)
!
! Compute a new block of the HH matrix
!
      IJ=(NDIMH+NDIMH**2)/2
      L1=NDIMH+1
      NDIMH=NDIMH+NTRIAL
      DO I=L1,NDIMH
       JST=1
       DO J=1,I
        IJ=IJ+1
        HH(IJ)=DDOT_(NDIM,HC(IST),1,C(JST),1)
        JST=JST+NDIM
       END DO
       IST=IST+NDIM
      END DO
!
      END SUBROUTINE HMAT
