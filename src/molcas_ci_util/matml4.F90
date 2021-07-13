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

subroutine MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL,ITRNSP)
! MULTIPLY A AND B TO GIVE C
!
! C = A * B             FOR ITRNSP = 0
! C = A(TRANSPOSED) * B FOR ITRNSP = 1
! C = A * B(TRANSPOSED) FOR ITRNSP = 2

implicit real*8(A-H,O-Z)
dimension A(NAROW,NACOL)
dimension B(NBROW,NBCOL)
dimension C(NCROW,NCCOL)

IZERO = 0
if ((NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL) == 0) IZERO = 1

if (ITRNSP == 0) then
  if (IZERO == 1) then
    call DCOPY_(NCROW*NCCOL,[0.0d0],0,C,1)
    do J=1,NCCOL
      do K=1,NBROW
        BKJ = B(K,J)
        call DAXPY_(NCROW,BKJ,A(1,K),1,C(1,J),1)
      end do
    end do
  else
    call DGEMM_('N','N',NCROW,NCCOL,NACOL,1.0d0,A,NAROW,B,NBROW,0.0d0,C,NCROW)
  end if
end if

if (ITRNSP == 1) then
  if (IZERO == 1) then
    do J=1,NCCOL
      do I=1,NCROW
        C(I,J) = DDOT_(NBROW,A(1,I),1,B(1,J),1)
      end do
    end do
  else
    call DGEMM_('T','N',NCROW,NCCOL,NAROW,1.0d0,A,NAROW,B,NBROW,0.0d0,C,NCROW)
  end if
end if

if (ITRNSP == 2) then
  if (IZERO == 1) then
    call DCOPY_(NCROW*NCCOL,[0.0d0],0,C,1)
    do J=1,NCCOL
      do K=1,NBCOL
        BJK = B(J,K)
        call DAXPY_(NCROW,BJK,A(1,K),1,C(1,J),1)
      end do
    end do
  else
    call DGEMM_('N','T',NCROW,NCCOL,NACOL,1.0d0,A,NAROW,B,NBROW,0.0d0,C,NCROW)
  end if
end if

return

end subroutine MATML4
