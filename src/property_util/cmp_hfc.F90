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
! Copyright (C) 2019, Thomas J. Duignan                                *
!               2021, Rulin Feng                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine dots the UHF spin density with the property integrals   *
! and trace the resulting matrix for each of the 9 components of the   *
! hyperfine magnetic integrals to obtain the 3 by 3 HFC tensor matrix  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
! nb     -  Number of total basis functions, input.                    *
! nat    -  Number of atoms, input.                                    *
!                                                                      *
!***********************************************************************

subroutine cmp_hfc(nb,nat)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nb, nat
integer(kind=iwp) :: iat, icomp, idir, iopt, ip, irc, jdir, jp, kdir, lu_one, nbtri, toper
real(kind=wp) :: amat(3,3), hfc(3,3), trace
character(len=8) :: label
real(kind=wp), allocatable :: sa(:,:), sd(:,:), sr(:,:), ta(:), td(:)

! Sizes of the matrices
nbtri = nb*(nb+1)/2
irc = -1
iopt = 0
lu_one = 2
toper = 255

call mma_allocate(td,nbtri,label='td')
call mma_allocate(sd,nb,nb,label='sd')
call get_dArray_chk('D1sao',td,nbtri)
call square(td,sd,nb,1,nb)

do ip=1,nb
  do jp=1,nb
    if (ip /= jp) then
      sd(jp,ip) = Half*sd(jp,ip)
    end if
  end do
end do

call mma_allocate(ta,nbtri+4,label='ta')
call mma_allocate(sa,nb,nb,label='sa')
call mma_allocate(sr,nb,nb,label='sr')

irc = -1
write(label,'(A,I3)') 'DEBUG',0
call opnone(irc,iopt,'ONEINT',lu_one)
if (irc /= 0) call Error()

do iat=1,nat
  icomp = 0
  do idir=1,3
    do jdir=1,3
      icomp = icomp+1
      trace = 0
      write(label,'(A,I3)') 'MAGXP',iat
      irc = -1
      call rdone(irc,iopt,label,icomp,ta,toper)
      if (irc /= 0) call Error()
      call square(ta,sa,nb,1,nb)
      call dgemm_('N','N',nb,nb,nb,One,sd,nb,sa,nb,Zero,sr,nb)
      do kdir=1,nb
        trace = trace+sr(kdir,kdir)
      end do
      hfc(idir,jdir) = trace
    end do
  end do

  ! spin field reduction to get the amat elements match
  amat(1,1) = hfc(2,2)+hfc(3,3)
  amat(2,2) = hfc(1,1)+hfc(3,3)
  amat(3,3) = hfc(1,1)+hfc(2,2)
  do idir=1,3
    do jdir=1,3
      if (idir /= jdir) amat(idir,jdir) = -hfc(jdir,idir)
    end do
  end do

  write(u6,*) ''
  write(u6,*) ''
  write(u6,'(A,I3)') 'Hyperfine coupling tensor matrix for atom:',iat
  write(u6,*) ''
  write(u6,'(A)') '   ---------------------------------------------------------'
  do idir=1,3
    write(u6,'(3ES20.10)') (-amat(idir,jdir),jdir=1,3)
  end do
  write(u6,'(A)') '   ---------------------------------------------------------'
end do
call Add_Info('AMAT',AMAT,9,5)

call mma_deallocate(td)
call mma_deallocate(sd)
call mma_deallocate(ta)
call mma_deallocate(sa)
call mma_deallocate(sr)
call clsone(irc,iopt)

return

contains

subroutine Error()
  write(u6,*) ' *** Error in subroutine cmp_hfc ***'
  write(u6,'(A,A)') '     Label = ',Label
  call Abend()
end subroutine Error

end subroutine cmp_hfc
