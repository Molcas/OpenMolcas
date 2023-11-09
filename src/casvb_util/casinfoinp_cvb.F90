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

subroutine casinfoinp_cvb()

use casvb_global, only: iorclos_d, iorcore_d, iorocc_d, istms2_d, istnel_d, istsy_d, mxirrep_ci, nstats_d, nstsym_d, weight_d
use Constants, only: Zero, One
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), parameter :: ncmp = 4, nstrin = 6
integer(kind=iwp) :: istr, nactel(3), nread
character(len=8), parameter :: string(nstrin) = ['RAS2    ','INACTIVE','FROZEN  ','NACTEL  ','SPIN    ','SYMMETRY']

do
  call fstring_cvb(string,nstrin,istr,ncmp,2)

  if (istr == 1) then
    ! 'RAS2'
    iorocc_d(:) = 0
    call int_cvb(iorocc_d,mxirrep_ci,nread,1)
  else if (istr == 2) then
    ! 'INACTIVE'
    iorclos_d(:) = 0
    call int_cvb(iorclos_d,mxirrep_ci,nread,1)
  else if (istr == 3) then
    ! 'FROZEN'
    iorcore_d(:) = 0
    call int_cvb(iorcore_d,mxirrep_ci,nread,1)
  end if
  if ((istr == 4) .or. (istr == 5) .or. (istr == 6)) then
    ! 'NACTEL' or 'SPIN' or 'SYMMETRY'
    if (nstsym_d == 0) then
      istnel_d(:) = 0
      istsy_d(:) = 0
      istms2_d(:) = 0
      nstats_d(:) = 0
      weight_d(:,:) = Zero
      nstsym_d = 1
      istsy_d(nstsym_d) = 1
      nstats_d(nstsym_d) = 1
      weight_d(1,nstsym_d) = One
    end if
  end if
  if (istr == 4) then
    ! 'NACTEL'
    nactel(:) = 0
    call int_cvb(nactel,3,nread,1)
    if ((nactel(2) /= 0) .or. (nactel(3) /= 0)) then
      write(u6,*) ' Illegal NACTEL read :',nactel
      write(u6,*) ' Use CASVB only with CASSCF wavefunctions!'
      call abend_cvb()
    end if
    istnel_d(nstsym_d) = nactel(1)
  else if (istr == 5) then
    ! 'SPIN'
    call int_cvb(istms2_d(nstsym_d),1,nread,1)
    istms2_d(nstsym_d) = istms2_d(nstsym_d)-1
  else if (istr == 6) then
    ! 'SYMMETRY'
    call int_cvb(istsy_d(nstsym_d),1,nread,1)
  end if
  ! Unrecognized keyword -- end of input:
  if (istr == 0) exit
end do

return

end subroutine casinfoinp_cvb
