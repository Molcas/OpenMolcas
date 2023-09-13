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
!************************************************
!** Subroutines to deal with CASSCF CI vectors **
!************************************************
!********************************
!** Routines involving CI only **
!********************************

subroutine cird_cvb(cvec,recn)
!***********************************************************************
!*                                                                     *
!*  CIRD   := Read CI vector.                                          *
!*                                                                     *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension cvec(*)
dimension idum(1)

ivec = nint(cvec(1))
iformat = iform_ci(ivec)
if (iformat == 0) then
  ioffs = 0
  call rdis_cvb(idum,1,recn,ioffs)
  iformat = idum(1)
  if (iformat /= iform_ci(ivec)) then
    write(6,*) ' Incompatible vector format on file.'
    write(6,*) ' Read :',iformat,' present :',iform_ci(ivec)
    call abend_cvb()
  end if
  call rdis_cvb(icnt_ci(ivec),1,recn,ioffs)
  call rdrs_cvb(work(iaddr_ci(ivec)),ndet,recn,ioffs)
else
  write(6,*) ' Unsupported format in CIRD :',iformat
  call abend_cvb()
end if

return

end subroutine cird_cvb
