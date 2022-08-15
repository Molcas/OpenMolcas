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
subroutine PDFTMemAlloc(mGrid,nOrbt)

use nq_pdft, only: d2RdRho2, d2RdRhodPi, d2ZdR2, dEdPi, dEdPiMO, dEdPix, dEdPiy, dEdPiz, dEdRho, dEdRhox, dEdRhoy, dEdRhoz, &
                   dF_dRhoamb, dF_dRhoapb, dF_dRhoxamb, dF_dRhoxapb, dF_dRhoyamb, dF_dRhoyapb, dF_dRhozamb, dF_dRhozapb, dRdPi, &
                   dRdRho, dRhodX, dRhodY, dRhodZ, dZdR, dZdRho, GdEdPiMO, GradPidFdRho, GradRdFdRho, GradRhodFdRho, MOas, MOax, &
                   MOay, MOaz, OneMz, OnePz, Pass1, Pass2, Pass3, RatioA, RhoAB, ZetaA
use KSDFT_Info, only: do_pdftpot
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid, nOrbt

call mma_allocate(RatioA,mGrid)
call mma_allocate(OnePz,mGrid)
call mma_allocate(OneMz,mGrid)
call mma_allocate(RhoAB,mGrid)
call mma_allocate(ZetaA,mGrid)
call mma_allocate(dZdR,mGrid)
call mma_allocate(Pass1,mGrid)
call mma_allocate(Pass2,mGrid)

! for ft-functional
call mma_allocate(Pass3,mGrid)

if (do_pdftPot) then
  call mma_allocate(dRdRho,mGrid)
  call mma_allocate(dRhodX,mGrid)
  call mma_allocate(dRhodY,mGrid)
  call mma_allocate(dRhodZ,mGrid)
  call mma_allocate(dF_dRhoapb,mGrid)
  call mma_allocate(dF_dRhoamb,mGrid)
  call mma_allocate(dF_dRhoxapb,mGrid)
  call mma_allocate(dF_dRhoyapb,mGrid)
  call mma_allocate(dF_dRhozapb,mGrid)
  call mma_allocate(dF_dRhoxamb,mGrid)
  call mma_allocate(dF_dRhoyamb,mGrid)
  call mma_allocate(dF_dRhozamb,mGrid)
  call mma_allocate(dEdRho,mGrid)
  call mma_allocate(dZdRho,mGrid)
  call mma_allocate(dEdRhox,mGrid)
  call mma_allocate(dEdRhoy,mGrid)
  call mma_allocate(dEdRhoz,mGrid)
  call mma_allocate(dEdPi,mGrid)
  call mma_allocate(GradRhodFdRho,mGrid)
  call mma_allocate(d2ZdR2,mGrid)
  call mma_allocate(d2RdRho2,mGrid)
  call mma_allocate(d2RdRhodPi,mGrid)
  call mma_allocate(MOas,mGrid,nOrbt)

  ! for ft-functional
  call mma_allocate(dRdPi,mGrid)
  call mma_allocate(GradRdFdRho,mGrid)
  call mma_allocate(GradPidFdRho,mGrid)
  call mma_allocate(dEdPix,mGrid)
  call mma_allocate(dEdPiy,mGrid)
  call mma_allocate(dEdPiz,mGrid)
  call mma_allocate(dEdPiMO,mGrid,nOrbt)
  call mma_allocate(GdEdPiMO,mGrid,nOrbt)
  call mma_allocate(MOax,mGrid,nOrbt)
  call mma_allocate(MOay,mGrid,nOrbt)
  call mma_allocate(MOaz,mGrid,nOrbt)
end if

end subroutine PDFTMemAlloc
