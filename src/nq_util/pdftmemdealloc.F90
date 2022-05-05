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
      Subroutine PDFTMemDeAlloc()

      use nq_pdft
      use KSDFT_Info, only: do_pdftpot
#include "stdalloc.fh"




      CALL mma_deallocate(RatioA        )
      CALL mma_deallocate(OnePz         )
      CALL mma_deallocate(OneMz         )
      CALL mma_deallocate(RhoAB         )
      CALL mma_deallocate(ZetaA         )
      CALL mma_deallocate(dZdR          )
      CALL mma_deallocate(Pass1         )
      CALL mma_deallocate(Pass2         )

!     for ft-functional
      CALL mma_deallocate(Pass3         )

      IF(do_pdftPot) THEN
       CALL mma_deallocate(dRdRho       )
       CALL mma_deallocate(dRhodX       )
       CALL mma_deallocate(dRhodY       )
       CALL mma_deallocate(dRhodZ       )
       CALL mma_deallocate(dF_dRhoapb   )
       CALL mma_deallocate(dF_dRhoamb   )
       CALL mma_deallocate(dF_dRhoxapb  )
       CALL mma_deallocate(dF_dRhoyapb  )
       CALL mma_deallocate(dF_dRhozapb  )
       CALL mma_deallocate(dF_dRhoxamb  )
       CALL mma_deallocate(dF_dRhoyamb  )
       CALL mma_deallocate(dF_dRhozamb  )
       CALL mma_deallocate(dEdRho       )
       CALL mma_deallocate(dZdRho       )
       CALL mma_deallocate(dEdRhox      )
       CALL mma_deallocate(dEdRhoy      )
       CALL mma_deallocate(dEdRhoz      )
       CALL mma_deallocate(dEdPi        )
       CALL mma_deallocate(GradRhodFdRho)
       CALL mma_deallocate(d2ZdR2       )
       CALL mma_deallocate(d2RdRho2     )
       CALL mma_deallocate(d2RdRhodPi   )
       CALL mma_deallocate(MOas         )
!      for ft-functional
       CALL mma_deallocate(dRdPi        )
       CALL mma_deallocate(GradRdFdRho  )
       CALL mma_deallocate(GradPidFdRho )
       CALL mma_deallocate(dEdPix       )
       CALL mma_deallocate(dEdPiy       )
       CALL mma_deallocate(dEdPiz       )
       CALL mma_deallocate(dEdPiMO      )
       CALL mma_deallocate(GdEdPiMO     )
       CALL mma_deallocate(MOax         )
       CALL mma_deallocate(MOay         )
       CALL mma_deallocate(MOaz         )
      END IF
      End Subroutine
