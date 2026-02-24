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
      SUBROUTINE INTCTL2(CMO,nCMO,DREF,nDREF,FIFA,NFIFA,HONE,nHONE)
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_grad, nStpGrd, FIMO_all, FIFA_all
      use caspt2_global, only: FIMO, FAMO
      use PrintLevel, only: DEBUG
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nBTri
      use definitions, only: iwp, wp
      IMPLICIT None
      integer(kind=iwp), intent(in):: nCMO, nDREF, NFIFA, nHONE
      Real(kind=wp), intent(in):: CMO(nCMO), DREF(nDREF), HONE(nHONE)
      Real(kind=wp), intent(out):: FIFA(NFIFA)

      LOGICAL(KIND=IWP), parameter:: IF_TRNSF=.False.
      Real(kind=wp), Allocatable:: FFAO(:), FIAO(:), FAAO(:)
      Real(kind=wp) tmp

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

      Call TraCho2(CMO,nCMO,DREF,nDREF,FFAO,FIAO,FAAO,IF_TRNSF)

       tmp=FFAO(1)
       tmp=sqrt(tmp)

      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL2 back from TRACHO2.'
        CALL XFLUSH(6)
      END IF
* All extra allocations inside tracho2 should now be gone.

* For gradient calculation, it is good to have FIAO and FAAO
      IF (do_grad.or.nStpGrd.eq.2) THEN
        !! FFAO has one-electron Hamiltonian

        FIMO_all(1:NBTri)=FFAO(:) + FIAO(:)
        FIFA_all(1:NBTri)=FIMO_all(:) + FAAO(:)

      END IF

* Transform them to MO basis:
      FIMO(:)=Zero
      FAMO(:)=Zero

c Compute FIMO, FAMO, ...  to workspace:
      Call FMat_Cho(CMO,SIZE(CMO),FIAO,FAAO,
     &              HONE,SIZE(HONE),FIMO,SIZE(FIMO),
     &                              FAMO,SIZE(FAMO),
     &                              FIFA,nFIFA)

      Call mma_deallocate(FFAO)
      Call mma_deallocate(FIAO)
      Call mma_deallocate(FAAO)

      END SUBROUTINE INTCTL2
