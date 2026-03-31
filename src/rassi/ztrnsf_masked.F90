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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine ZTRNSF_MASKED(N,UR,UI,AR,AI,IJ,IST,INUM,JST,JNUM)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N, IJ(4), INUM, IST(INUM), JNUM, JST(JNUM)
real(kind=wp) :: UR(N,N), UI(N,N), AR(N,N), AI(N,N)
integer(kind=iwp) :: I, II, J, JJ, NI, NJ
real(kind=wp), allocatable :: MI(:,:), MR(:,:), TI(:,:), TR(:,:), VI(:,:), VR(:,:)

NI = IJ(2)-IJ(1)+1
NJ = IJ(4)-IJ(3)+1

call mma_allocate(MR,INUM,JNUM,LABEL='MR')
call mma_allocate(MI,INUM,JNUM,LABEL='MI')
do J=1,JNUM
  do I=1,INUM
    MR(I,J) = AR(IST(I),JST(J))
    MI(I,J) = AI(IST(I),JST(J))
  end do
end do

call mma_allocate(VR,JNUM,NJ,LABEL='VR')
call mma_allocate(VI,JNUM,NJ,LABEL='VI')
do J=1,NJ
  JJ = IJ(3)+J-1
  do I=1,JNUM
    VR(I,J) = UR(JST(I),JJ)
    VI(I,J) = UI(JST(I),JJ)
  end do
end do

call mma_allocate(TR,INUM,NJ,LABEL='TR')
call mma_allocate(TI,INUM,NJ,LABEL='TI')
call DGEMM_('N','N',INUM,NJ,JNUM,One,MR,INUM,VR,JNUM,Zero,TR,INUM)
call DGEMM_('N','N',INUM,NJ,JNUM,-One,MI,INUM,VI,JNUM,One,TR,INUM)
call DGEMM_('N','N',INUM,NJ,JNUM,One,MR,INUM,VI,JNUM,Zero,TI,INUM)
call DGEMM_('N','N',INUM,NJ,JNUM,One,MI,INUM,VR,JNUM,One,TI,INUM)

call mma_deallocate(VR)
call mma_deallocate(VI)
call mma_allocate(VR,INUM,NI,LABEL='VR')
call mma_allocate(VI,INUM,NI,LABEL='VI')
do J=1,NI
  JJ = IJ(1)+J-1
  do I=1,INUM
    VR(I,J) = UR(IST(I),JJ)
    VI(I,J) = UI(IST(I),JJ)
  end do
end do

call mma_deallocate(MR)
call mma_deallocate(MI)
call mma_allocate(MR,NI,NJ,LABEL='MR')
call mma_allocate(MI,NI,NJ,LABEL='MI')
call DGEMM_('T','N',NI,NJ,INUM,One,VR,INUM,TR,INUM,Zero,MR,NI)
call DGEMM_('T','N',NI,NJ,INUM,One,VI,INUM,TI,INUM,One,MR,NI)
call DGEMM_('T','N',NI,NJ,INUM,One,VR,INUM,TI,INUM,Zero,MI,NI)
call DGEMM_('T','N',NI,NJ,INUM,-One,VI,INUM,TR,INUM,One,MI,NI)

call DCOPY_(N*N,[Zero],0,AR,1)
call DCOPY_(N*N,[Zero],0,AI,1)
do J=1,NJ
  JJ = IJ(3)+J-1
  do I=1,NI
    II = IJ(1)+I-1
    AR(II,JJ) = MR(I,J)
    AI(II,JJ) = MI(I,J)
  end do
end do

call mma_deallocate(TR)
call mma_deallocate(TI)
call mma_deallocate(VR)
call mma_deallocate(VI)
call mma_deallocate(MR)
call mma_deallocate(MI)

return

end subroutine ZTRNSF_MASKED
