************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE INTCTL2(IF_TRNSF)
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_grad, nStpGrd, FIMO_all, FIFA_all
      use caspt2_global, only: CMO, FIMO, FAMO, HONE, DREF
      use PrintLevel, only: DEBUG
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nBTri
      IMPLICIT None
      LOGICAL IF_TRNSF

      Real*8, Allocatable:: FFAO(:), FIAO(:), FAAO(:)

* Compute using Cholesky vectors.
* Frozen, inactive and active Fock matrix in AO basis:
      Call mma_allocate(FFAO,NBTRI,LABEL='FFAO')
      Call mma_allocate(FIAO,NBTRI,LABEL='FIAO')
      Call mma_allocate(FAAO,NBTRI,LABEL='FAAO')
* tracho2 makes many allocations but should deallocate everything
* before its return.
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL2 calling TRACHO2...'
        CALL XFLUSH(6)
      END IF
      Call TraCho2(CMO,SIZE(CMO),DREF,SIZE(DREF),FFAO,FIAO,FAAO,
     &             IF_TRNSF)
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL2 back from TRACHO2.'
        CALL XFLUSH(6)
      END IF
* All extra allocations inside tracho2 should now be gone.

* For gradient calculation, it is good to have FIAO and FAAO
      IF (do_grad.or.nStpGrd.eq.2) THEN
        !! FFAO has one-electron Hamiltonian
        CALL DCOPY_(NBTRI,FFAO,1,FIMO_all,1)
        CALL DAXPY_(NBTRI,One,FIAO,1,FIMO_all,1)
        CALL DCOPY_(NBTRI,FIMO_all,1,FIFA_all,1)
        CALL DAXPY_(NBTRI,One,FAAO,1,FIFA_all,1)
      END IF
* Transform them to MO basis:
      HONE(:)=Zero
      FIMO(:)=Zero
      FAMO(:)=Zero
c Compute FIMO, FAMO, ...  to workspace:
      Call FMat_Cho(CMO,SIZE(CMO),FFAO,FIAO,FAAO,
     &              HONE,SIZE(HONE),FIMO,SIZE(FIMO),FAMO,SIZE(FAMO))

      Call mma_deallocate(FFAO)
      Call mma_deallocate(FIAO)
      Call mma_deallocate(FAAO)

      END SUBROUTINE INTCTL2
