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

function recinpcmp_cvb(ifield)

implicit real*8(a-h,o-z)
logical recinpcmp_cvb
! ... Files/Hamiltonian available ...
logical, external :: valid_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
logical done

if (.not. valid_cvb(recinp_old)) then
  recinpcmp_cvb = .true.
else
  call rdioff_cvb(ifield,recinp,ioff1)
  call rdioff_cvb(ifield+1,recinp,ioff2)
  call rdioff_cvb(ifield,recinp_old,joff1)
  call rdioff_cvb(ifield+1,recinp_old,joff2)
  if (ioff2-ioff1 /= joff2-joff1) then
    recinpcmp_cvb = .true.
  else
    i1 = mstackr_cvb(ioff2-ioff1)
    j1 = mstackr_cvb(joff2-joff1)
    call rdr_cvb(work(i1),ioff2-ioff1,recinp,ioff1)
    call rdr_cvb(work(j1),joff2-joff1,recinp_old,joff1)
    done = .false.
    do i=0,ioff2-ioff1-1
      if (work(i+i1) /= work(i+j1)) then
        recinpcmp_cvb = .true.
        done = .true.
        exit
      end if
    end do
    if (.not. done) recinpcmp_cvb = .false.
    call mfreer_cvb(i1)
  end if
end if

return

end function recinpcmp_cvb
