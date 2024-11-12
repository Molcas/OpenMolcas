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
      SUBROUTINE ZTRNSF(N,UR,UI,AR,AI)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer N
      Real*8 UR(N,N),UI(N,N)
      Real*8 AR(N,N),AI(N,N)

      Real*8, Allocatable:: CR(:,:), CI(:,:)

      Call mma_allocate(CR,N,N,Label='CR')
      Call mma_allocate(CI,N,N,Label='CI')
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AR,N,UR,N,0.0D0,CR,N)
      CALL DGEMM_('N','N',N,N,N,-1.0D0,AI,N,UI,N,1.0D0,CR,N)
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AR,N,UI,N,0.0D0,CI,N)
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AI,N,UR,N,1.0D0,CI,N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UR,N,CR,N,0.0D0,AR,N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UI,N,CI,N,1.0D0,AR,N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UR,N,CI,N,0.0D0,AI,N)
      CALL DGEMM_('T','N',N,N,N,-1.0D0,UI,N,CR,N,1.0D0,AI,N)
      Call mma_deallocate(CI)
      Call mma_deallocate(CR)

      END SUBROUTINE ZTRNSF
