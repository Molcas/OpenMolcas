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

subroutine symelminp_cvb(rsymelm,nsyme,tags,izeta,mxirrep,mxorb,mxsyme,ityp)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nsyme, mxirrep, mxorb, mxsyme, ityp(mxorb)
real(kind=wp), intent(inout) :: rsymelm(mxorb,mxorb,nsyme)
character(len=3), intent(inout) :: tags(mxsyme)
integer(kind=iwp), intent(inout) :: izeta(*)
integer(kind=iwp) :: i, i_dim, iaux(1), io, iorb, irrep, isgn, istr2, jor, jorb, nread
real(kind=wp) :: daux(1)
integer(kind=iwp), allocatable :: tmp(:)
integer(kind=iwp), parameter :: ncmp = 4, nsign = 2, nsymelm = 5
character(len=*), parameter :: sgn(nsign) = ['+       ','-       '], &
                               symelm(nsymelm) = ['IRREPS  ','COEFFS  ','TRANS   ','END     ','ENDSYMEL']
logical(kind=iwp), external :: mxorth_cvb ! ... Matrices (orthogonal/determinant) ...

tags(nsyme) = ' '
call string_cvb(tags(nsyme),1,nread,1)
call fstring_cvb(sgn,nsign,isgn,ncmp,1)
if (isgn == 1) then
  izeta(nsyme) = 1
else if (isgn == 2) then
  izeta(nsyme) = -1
else
  izeta(nsyme) = 0
end if
call unitmat(rsymelm(:,:,nsyme),mxorb)

do
  call fstring_cvb(symelm,nsymelm,istr2,ncmp,2)
  if (istr2 == 1) then
    ! 'IRREPS'
    do i=1,mxirrep
      iaux = 0
      call int_cvb(iaux,1,nread,0)
      irrep = iaux(1)
      if (irrep /= 0) then
        do iorb=1,mxorb
          if (irrep == ityp(iorb)) rsymelm(iorb,iorb,nsyme) = -One
        end do
      end if
    end do
  else if (istr2 == 2) then
    ! 'COEFFS'
    do i=1,mxorb
      iaux = 0
      call int_cvb(iaux,1,nread,0)
      iorb = iaux(1)
      if (iorb /= 0) then
        rsymelm(iorb,iorb,nsyme) = -One
      else
        exit
      end if
    end do
  else if (istr2 == 3) then
    ! 'TRANS'
    iaux = 0
    call int_cvb(iaux,1,nread,0)
    i_dim = iaux(1)
    if ((i_dim < 1) .or. (i_dim > mxorb)) then
      write(u6,*) ' Illegal dimension in TRANS:',i_dim,mxorb
      call abend_cvb()
    end if
    call mma_allocate(tmp,i_dim,label='tmp')
    do i=1,i_dim
      call int_cvb(iaux,1,nread,0)
      iorb = iaux(1)
      if ((iorb < 1) .or. (iorb > mxorb)) then
        write(u6,*) ' Illegal orbital number in TRANS:',iorb
        call abend_cvb()
      end if
      tmp(i) = iorb
    end do
    do io=1,i_dim
      iorb = tmp(io)
      do jor=1,i_dim
        jorb = tmp(jor)
        daux = Zero
        call real_cvb(daux(1),1,nread,0)
        rsymelm(iorb,jorb,nsyme) = daux(1)
      end do
    end do
    call mma_deallocate(tmp)
  end if
  ! 'END' , 'ENDSYMEL' or unrecognized keyword -- end SYMELM input:
  if ((istr2 == 4) .or. (istr2 == 5) .or. (istr2 == 0)) exit
end do
if (.not. mxorth_cvb(rsymelm(:,:,nsyme),mxorb)) then
  write(u6,*) ' Symmetry element ',tags(nsyme),' not orthogonal!'
  write(u6,*) ' Check usage of TRANS keyword.'
  call abend_cvb()
end if

return

end subroutine symelminp_cvb
