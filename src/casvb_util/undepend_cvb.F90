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

subroutine undepend_cvb(chr1,chr2)

use casvb_global, only: charobj, i_dep_on_j, ioffs, iprint, j_dep_on_i, joffs, mustdeclare, ndep_ij, ndep_ji, nobj

implicit real*8(a-h,o-z)
character*(*) chr1, chr2
logical done

ic = 3

do
  iobj = 0
  jobj = 0
  do i=1,nobj
    if (charobj(i) == chr1) iobj = i
    if (charobj(i) == chr2) jobj = i
  end do
  if (iobj == 0) then
    if (mustdeclare) then
      write(6,*) ' Make object not found :',chr1
      call abend_cvb()
    end if
    call decl_cvb(chr1)
  else if (jobj == 0) then
    if (mustdeclare) then
      write(6,*) ' Make object not found :',chr2
      call abend_cvb()
    end if
    call decl_cvb(chr2)
  else
    exit
  end if
end do

if (iprint >= 10) write(6,*) ' Cancel I depends on J :',iobj,jobj
n_cancelled = 0
if (mod(ic,2) == 1) then
  do
    done = .false.
    do i=ioffs(iobj)+1,ioffs(iobj+1)
      if (i_dep_on_j(i) == jobj) then
        do j=i,ioffs(nobj+1)-1
          i_dep_on_j(j) = i_dep_on_j(j+1)
        end do
        do ii=iobj+1,nobj+1
          ioffs(ii) = ioffs(ii)-1
        end do
        n_cancelled = n_cancelled+1
        done = .true.
        exit
      end if
    end do
    if (.not. done) exit
  end do
end if

m_cancelled = 0
if (ic >= 2) then
  do
    done = .false.
    do i=joffs(jobj)+1,joffs(jobj+1)
      if (j_dep_on_i(i) == iobj) then
        do j=i,joffs(nobj+1)-1
          j_dep_on_i(j) = j_dep_on_i(j+1)
        end do
        do ii=jobj+1,nobj+1
          joffs(ii) = joffs(ii)-1
        end do
        m_cancelled = m_cancelled+1
        done = .true.
        exit
      end if
    end do
    if (.not. done) exit
  end do
end if

ndep_ij = ndep_ij-n_cancelled
ndep_ji = ndep_ji-m_cancelled

return

end subroutine undepend_cvb
