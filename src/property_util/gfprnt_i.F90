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

subroutine GFPrnt_i(EVal,nDim)

implicit real*8(a-h,o-z)
real*8 EVal(nDim)
character*80 format, Line*120

LUt = 6
Inc = 6
do iHarm=1,nDim,Inc
  Jnc = min(Inc,nDim-iHarm+1)
  write(format,'(A,I3,A)') '(5X,A10,1x,',Jnc,'I10)'
  write(LUt,format) ' ',(i,i=iHarm,iHarm+Jnc-1)
  write(LUt,*)

  write(format,'(A,I3,A)') '(A12,1x,',Jnc,'F10.2)'
  Line = ' '
  write(Line,format) 'Freq.',(EVal(i),i=iHarm,iHarm+Jnc-1)
  do i=1,120
    if (Line(i:i) == '-') Line(i:i) = 'i'
  end do
  write(LUt,'(5X,A)') Line
  write(LUt,*)

  write(LUt,*)
end do

return

end subroutine GFPrnt_i
