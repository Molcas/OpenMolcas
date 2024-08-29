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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

module mspdftgrad
  use definitions,only:wp
  implicit none
  private

  real(kind=wp),allocatable :: F1MS(:,:),F2MS(:,:),FocMS(:,:),FxyMS(:,:),P2Mot(:,:)

  public :: F1MS,F2MS,FocMS,FxyMS,P2Mot

  public :: mspdftgrad_init,mspdftgrad_free

contains

  subroutine mspdftgrad_init()
    use constants,only:zero
    use stdalloc,only:mma_allocate
    use rasscf_global,only:nroots,nacpr2,nTot4
    implicit none

#include "rasdim.fh"
#include "general.fh"

    call mma_allocate(F1MS,nTot1,nRoots,label="F1MS")
    call mma_allocate(F2MS,NacPR2,nRoots,label="F2MS")
    call mma_allocate(FocMS,nTot1,nroots,label="FocMS")
    call mma_allocate(FxyMS,nTot4,nRoots,label="FxyMS")
    call mma_allocate(P2MOT,nacpr2,nroots,label="P2MO")

    P2MOT = zero

  endsubroutine

  subroutine mspdftgrad_free()
    use stdalloc,only:mma_deallocate
    implicit none

    ! Should in theory check before deallocating!
    call mma_deallocate(F1MS)
    call mma_deallocate(F2MS)
    call mma_deallocate(FocMS)
    call mma_deallocate(FxyMS)
    call mma_deallocate(P2MOT)
  endsubroutine

endmodule
