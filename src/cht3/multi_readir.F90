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

subroutine multi_readir(G,lg,ifile,ias)
! Direct fortran I/O with irregular data records
!
! Each data record is assumed to begin aligned with a disk record
! but may span several disk records.
!
! The price to pay for this flexibility is
! 1) remembering the correspondance between data record numbers
!    and corresponding initial disk records.
! 2) wasting some disk space, in average about nu*nblock/2
!    where nu is the number of data records and nblock the disk
!    record size.
!
! Arguments
!    G        Buffer (real*8 words)
!    lg       Buffer length
!    ifile    file unit
!    ias      direct access record to start with
!    (nblock  direct access record length, defined in module)
!
! PV/LAOG, 22 may 2003.

use ChT3_global, only: IOPT, nblock
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lg, ifile, ias
real(kind=wp), intent(out) :: G(lg)
integer(kind=iwp) :: iloc, irest, k, kas

iloc = 1
irest = lg
kas = ias

do while (irest > 0)
  k = min(irest,nblock)
  if (kas <= iopt(2)) then
    read(ifile,rec=kas) G(iloc:iloc+k-1)
  else
    read(ifile+1,rec=kas-iopt(2)) G(iloc:iloc+k-1)
  end if
  iloc = iloc+k
  irest = irest-k
  kas = kas+1
end do

end subroutine multi_readir
