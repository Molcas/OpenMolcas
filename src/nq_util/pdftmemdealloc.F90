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
subroutine PDFTMemDeAlloc()

use nq_pdft, only: d2RdRho2, d2RdRhodPi, d2ZdR2, dEdPi, dEdPiMO, dEdPix, dEdPiy, dEdPiz, dEdRho, dEdRhox, dEdRhoy, dEdRhoz, &
                   dF_dRhoamb, dF_dRhoapb, dF_dRhoxamb, dF_dRhoxapb, dF_dRhoyamb, dF_dRhoyapb, dF_dRhozamb, dF_dRhozapb, dRdPi, &
                   dRdRho, dRhodX, dRhodY, dRhodZ, dZdR, dZdRho, GdEdPiMO, GradPidFdRho, GradRdFdRho, GradRhodFdRho, MOas, MOax, &
                   MOay, MOaz, OneMz, OnePz, Pass1, Pass2, Pass3, RatioA, RhoAB, ZetaA
use KSDFT_Info, only: do_pdftpot
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(RatioA)
call mma_deallocate(OnePz)
call mma_deallocate(OneMz)
call mma_deallocate(RhoAB)
call mma_deallocate(ZetaA)
call mma_deallocate(dZdR)
call mma_deallocate(Pass1)
call mma_deallocate(Pass2)

! for ft-functional
call mma_deallocate(Pass3)

if (do_pdftPot) then
  call mma_deallocate(dRdRho)
  call mma_deallocate(dRhodX)
  call mma_deallocate(dRhodY)
  call mma_deallocate(dRhodZ)
  call mma_deallocate(dF_dRhoapb)
  call mma_deallocate(dF_dRhoamb)
  call mma_deallocate(dF_dRhoxapb)
  call mma_deallocate(dF_dRhoyapb)
  call mma_deallocate(dF_dRhozapb)
  call mma_deallocate(dF_dRhoxamb)
  call mma_deallocate(dF_dRhoyamb)
  call mma_deallocate(dF_dRhozamb)
  call mma_deallocate(dEdRho)
  call mma_deallocate(dZdRho)
  call mma_deallocate(dEdRhox)
  call mma_deallocate(dEdRhoy)
  call mma_deallocate(dEdRhoz)
  call mma_deallocate(dEdPi)
  call mma_deallocate(GradRhodFdRho)
  call mma_deallocate(d2ZdR2)
  call mma_deallocate(d2RdRho2)
  call mma_deallocate(d2RdRhodPi)
  call mma_deallocate(MOas)
  ! for ft-functional
  call mma_deallocate(dRdPi)
  call mma_deallocate(GradRdFdRho)
  call mma_deallocate(GradPidFdRho)
  call mma_deallocate(dEdPix)
  call mma_deallocate(dEdPiy)
  call mma_deallocate(dEdPiz)
  call mma_deallocate(dEdPiMO)
  call mma_deallocate(GdEdPiMO)
  call mma_deallocate(MOax)
  call mma_deallocate(MOay)
  call mma_deallocate(MOaz)
end if

end subroutine PDFTMemDeAlloc
