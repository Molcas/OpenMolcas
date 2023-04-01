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

subroutine multc0(wrk,wrksize,mvec,ix,mapdc,key)
! This routine realizes multiplying according mvec
! for C=A*B
! N.B. if key=0, C file is not vanished (i.e. can be used for
! adding to some existing file)
!
! If C=A*B process is faster or comparable with C=AT*B then mchntyp should be set to 1.
! If C=AT*B is significantly faster than C=A*T (more than 20%), than mchntyp should be set to 2. (default is 1)
! if mchntyp is 2, then
! 1) processes with scale(A)/scale(B) > scalelim will be calculated as C=A*B
! 2) processes with scale(A)/scale(B) < scalelim will be calculated as C=AT*B
! Note, that for mchntyp =2 more memory is required, due to requirement of
! aditional o2v2 help file possd0 (parameter possd0 is transported through ccsd1.fh, not
! through ccsd2.fh)

#include "ccsd1.fh"
#include "wrk.fh"
integer mvec(1:4096,1:7)
integer ix, key
integer mapdc(0:512,1:6)
! help variables
integer nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6
integer iix, ic
real*8 scale_

!1 set C=0

if (key == 1) then

  ! C matrix must be vanished

  do ic=1,mapdc(0,5)
    nhelp1 = mapdc(ic,1)
    nhelp2 = mapdc(ic,2)
    call mv0zero(nhelp2,nhelp2,wrk(nhelp1))
  end do

end if

!2 C=C+A*B

if (ix == 0) then
  return
end if

do iix=1,ix

  ! skip this summation if yes/no=0
  if (mvec(iix,1) == 0) goto 200

  ! realize individial summation

  ! def positions of A,B,C
  nhelp1 = mvec(iix,2)
  nhelp2 = mvec(iix,3)
  nhelp3 = mvec(iix,4)

  ! def rowA(rowC), colA(rowB,sum), colB(colC)
  nhelp4 = mvec(iix,5)
  nhelp5 = mvec(iix,6)
  nhelp6 = mvec(iix,7)

  if (mchntyp == 1) then

    !* Typ 1
    call mc0c1a3b(nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

  else

    !* Typ2
    scale_ = (1.0d0*nhelp4)/(1.0d0*nhelp6)
    if (scale_ > slim) then
      call mc0c1a3b(nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

    else
      ! map D=AT
      call map21(wrk(nhelp1),wrk(possd0),nhelp4,nhelp5,2,1,1)
      ! calc C=DT*B
      call mc0c1at3b(nhelp5,nhelp4,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(possd0),wrk(nhelp2),wrk(nhelp3))
    end if

  end if

  200 continue
end do

return

end subroutine multc0
