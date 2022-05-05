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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************
      Subroutine PDFTMemAlloc(mGrid,nOrbt)

      use nq_pdft
      use KSDFT_Info, only: do_pdftpot
#include "stdalloc.fh"

      INTEGER mGrid,nOrbt


      CALL mma_allocate(RatioA        ,mGrid)
      CALL mma_allocate(OnePz         ,mGrid)
      CALL mma_allocate(OneMz         ,mGrid)
      CALL mma_allocate(RhoAB         ,mGrid)
      CALL mma_allocate(ZetaA         ,mGrid)
      CALL mma_allocate(dZdR          ,mGrid)
      CALL mma_allocate(Pass1         ,mGrid)
      CALL mma_allocate(Pass2         ,mGrid)

!     for ft-functional
      CALL mma_allocate(Pass3         ,mGrid)

      IF(do_pdftPot) THEN
       CALL mma_allocate(dRdRho       ,mGrid)
       CALL mma_allocate(dRhodX       ,mGrid)
       CALL mma_allocate(dRhodY       ,mGrid)
       CALL mma_allocate(dRhodZ       ,mGrid)
       CALL mma_allocate(dF_dRhoapb   ,mGrid)
       CALL mma_allocate(dF_dRhoamb   ,mGrid)
       CALL mma_allocate(dF_dRhoxapb  ,mGrid)
       CALL mma_allocate(dF_dRhoyapb  ,mGrid)
       CALL mma_allocate(dF_dRhozapb  ,mGrid)
       CALL mma_allocate(dF_dRhoxamb  ,mGrid)
       CALL mma_allocate(dF_dRhoyamb  ,mGrid)
       CALL mma_allocate(dF_dRhozamb  ,mGrid)
       CALL mma_allocate(dEdRho       ,mGrid)
       CALL mma_allocate(dZdRho       ,mGrid)
       CALL mma_allocate(dEdRhox      ,mGrid)
       CALL mma_allocate(dEdRhoy      ,mGrid)
       CALL mma_allocate(dEdRhoz      ,mGrid)
       CALL mma_allocate(dEdPi        ,mGrid)
       CALL mma_allocate(GradRhodFdRho,mGrid)
       CALL mma_allocate(d2ZdR2       ,mGrid)
       CALL mma_allocate(d2RdRho2     ,mGrid)
       CALL mma_allocate(d2RdRhodPi   ,mGrid)
       CALL mma_allocate(MOas         ,mGrid*nOrbt)

!      for ft-functional
       CALL mma_allocate(dRdPi        ,mGrid)
       CALL mma_allocate(GradRdFdRho  ,mGrid)
       CALL mma_allocate(GradPidFdRho ,mGrid)
       CALL mma_allocate(dEdPix       ,mGrid)
       CALL mma_allocate(dEdPiy       ,mGrid)
       CALL mma_allocate(dEdPiz       ,mGrid)
       CALL mma_allocate(dEdPiMO      ,mGrid*nOrbt)
       CALL mma_allocate(GdEdPiMO     ,mGrid*nOrbt)
       CALL mma_allocate(MOax         ,mGrid*nOrbt)
       CALL mma_allocate(MOay         ,mGrid*nOrbt)
       CALL mma_allocate(MOaz         ,mGrid*nOrbt)
      END IF
      End Subroutine
