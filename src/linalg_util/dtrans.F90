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

subroutine DTRANS(NROWS,NCOLS,A,LDA,B,LDB)
! Copy-transpose matrix A to B using blocks to optimize cache
! usage or extended BLAS functionality if present.
! B := transpose(A)
!
! double precision version

#include "intent.fh"

#ifdef _MKL_
use Constants, only: One
#endif
use Definitions, only: wp, iwp, u6

implicit none
! arguments
integer(kind=iwp), intent(in) :: NROWS, NCOLS, LDA, LDB
real(kind=wp), intent(in) :: A(LDA,*)
real(kind=wp), intent(_OUT_) :: B(LDB,*)
! local variables
#ifndef _MKL_
integer(kind=iwp) :: I, IB, J, JB, LCOLS, LROWS, MAXCOL, MAXROW, NBLKSZ
#endif

if ((NROWS <= 0) .or. (NCOLS <= 0)) then
  write(u6,'(1X,A)') 'DTRANS: Error: invalid dimension(s)'
  write(u6,'(1X,2(A,I9))') 'NROWS = ',NROWS,'NCOLS = ',NCOLS
  call AbEnd
else if ((NROWS > LDA) .or. (NCOLS > LDB)) then
  write(u6,'(1X,A)') 'DTRANS: Error: dimension(s) out-of-bounds'
  write(u6,'(1X,2(A,I9))') 'NROWS = ',NROWS,'NCOLS = ',NCOLS
  write(u6,'(1X,2(A,I9))') 'LDA   = ',LDA,'LDB   = ',LDB
  call AbEnd
end if

#ifdef _MKL_
call mkl_domatcopy('C','T',NROWS,NCOLS,One,A,LDA,B,LDB)
#else
NBLKSZ = 8
LROWS = mod(NROWS,NBLKSZ)
LCOLS = mod(NCOLS,NBLKSZ)
MAXROW = NROWS-LROWS
MAXCOL = NCOLS-LCOLS
if ((MAXROW > 0) .and. (MAXCOL > 0)) then
  do IB=1,MAXROW,NBLKSZ
    do JB=1,MAXCOL,NBLKSZ
      do I=IB,IB+NBLKSZ-1
        do J=JB,JB+NBLKSZ-1
          B(J,I) = A(I,J)
        end do
      end do
    end do
  end do
end if
! remainder of the blocks
if ((MAXROW > 0) .and. (LCOLS > 0)) then
  do IB=1,MAXROW,NBLKSZ
    do I=IB,IB+NBLKSZ-1
      do J=MAXCOL+1,MAXCOL+LCOLS
        B(J,I) = A(I,J)
      end do
    end do
  end do
end if
if ((MAXCOL > 0) .and. (LROWS > 0)) then
  do JB=1,MAXCOL,NBLKSZ
    do I=MAXROW+1,MAXROW+LROWS
      do J=JB,JB+NBLKSZ-1
        B(J,I) = A(I,J)
      end do
    end do
  end do
end if
! remainder block
if ((LROWS > 0) .and. (LCOLS > 0)) then
  do I=MAXROW+1,MAXROW+LROWS
    do J=MAXCOL+1,MAXCOL+LCOLS
      B(J,I) = A(I,J)
    end do
  end do
end if
#endif

end subroutine DTRANS
