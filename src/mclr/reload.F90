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
      Subroutine ReLoad(A,idsym,NBAS1,NBAS2)
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: ipMat, nDens2
      use input_mclr, only: nSym
      Implicit None
      Real*8 A(*)
      Integer idSym
      Integer nbas2(nsym),nbas1(nsym)
      Real*8, Allocatable:: ATemp(:)
      Integer iS, jS, j

      Call mma_allocate(ATemp,ndens2,Label='ATemp')

      Do iS=1,nsym
       js=ieor(is-1,idsym-1) +1
       if (min(nbas1(is),nbas2(is)) < 1) cycle
       Do j=0,Min(nbas2(js),nbas1(js))-1
        call dcopy_(Min(nbas1(is),nbas2(is)),                           &
     &             A(ipMat(is,js)+j*nbas1(is)),1,                       &
     &         ATemp(ipmat(is,js)+j*nbas2(is)),1)
       End Do
      End Do
      call dcopy_(ndens2,ATemp,1,A,1)
      Call mma_deallocate(ATemp)
      End Subroutine ReLoad
