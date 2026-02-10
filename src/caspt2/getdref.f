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
      SUBROUTINE GETDREF(DREF,NDREF)
      use definitions, onlY: iwp, wp, u6
      use constants, only: Zero
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT
      use pt2_guga, only: NG1
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: NDREF
      real(kind=wp), intent(out):: DREF(NDREF)

      real(kind=wp), ALLOCATABLE:: G1(:)
      integer(kind=iwp) I, J, IJ

* Get active 1-el density matrix GAMMA1 and
* construct DREF in a tringular storage.

* Remember: NDREF=1 if NASHT=0.
      DREF(1)=Zero
      IF(NASHT==0) RETURN
* Active 1-el density matrix:
      CALL mma_allocate(G1,NG1,Label='G1')
      CALL PT2_GET(NG1,'GAMMA1',G1)
      DO I=1,NASHT
        DO J=1,I
          IJ=(I*(I-1))/2+J
          DREF(IJ)=G1(I+NASHT*(J-1))
        END DO
      END DO
      CALL mma_deallocate(G1)

      IF(IPRGLB.GE.DEBUG) THEN
        WRITE(u6,*)' GETDREF has constructed DREF.'
      END IF

      END SUBROUTINE GETDREF

