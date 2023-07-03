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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine List(Line,Lbl,gq,mInt,nIter)
!***********************************************************************
!                                                                      *
! Object: to print gradient or internal coordinate lists               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             1993                                                     *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 gq(mInt,nIter)
character Lbl(mInt)*8, format*72, Line*(*)
#include "print.fh"
#include "real.fh"

Lu = 6

write(Lu,*)
write(Lu,*)
write(Lu,*) Line

MxWdth = 132
nLbl = 8+1
nRow = 10
inc = min((MxWdth-nLbl)/nRow,nIter)

do ii=1,nIter,inc
  write(Lu,*)
  write(format,'(A,I2,A)') '(A,1X,',inc,'(I5,5X))'
  write(Lu,format) 'Iter.   ',(i,i=ii,min(ii+inc-1,nIter))
  write(Lu,*)
  write(format,'(A,I2,A)') '(A,1X,',inc,'(F9.5,1X))'
  do igq=1,mInt
    write(Lu,format) Lbl(igq),(gq(igq,i),i=ii,min(ii+inc-1,nIter))
  end do
  write(Lu,*)
  write(Lu,*)
end do
write(Lu,*)

return

end subroutine List
