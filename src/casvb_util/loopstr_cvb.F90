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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine loopstr_cvb(iocc,indx,nel,norb)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nel, norb
integer(kind=iwp), intent(inout) :: iocc(nel), indx
integer(kind=iwp) :: iel, jel
logical(kind=iwp) :: done

indx = indx+1
! Find electron for which orbital number can be increased:
done = .false.
do iel=1,nel-1
  if (iocc(iel+1) > iocc(iel)+1) then
    done = .true.
    exit
  end if
end do
if (.not. done) then
  iel = nel
  if (iocc(iel) >= norb) then
    call loopstr0_cvb(iocc,indx,nel,norb)
    return
  end if
end if

iocc(iel) = iocc(iel)+1
do jel=1,iel-1
  iocc(jel) = jel
end do

return

end subroutine loopstr_cvb
