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

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nb, nat
integer(kind=iwp) :: iat, icomp, idir, iopt, ip, irc, isa, isd, isr, ita, itd, jdir, jp, kdir, lu_one, nb2, nbtri, stri, toper
real(kind=wp) :: amat(3,3), hfc(3,3), trace
character(len=8) :: label
#include "WrkSpc.fh"

! Sizes of the matrices
nb2 = nb*nb
nbtri = nb*(nb+1)/2
irc = -1
iopt = 0
lu_one = 2
toper = 255

call allocate_work(itd,nbtri+4)
call allocate_work(isd,nb2+4)
call get_d1sao(Work(itd),nbtri)
call square(Work(itd),Work(isd),nb,1,nb)

stri = 0
do ip=0,nb-1
  stri = isd+ip*nb
  do jp=0,nb-1
    if (ip /= jp) then
      Work(stri+jp) = Half*Work(stri+jp)
    end if
  end do
end do

call allocate_work(ita,nbtri+4)
call allocate_work(isa,nb2+4)
call allocate_work(isr,nb2+4)

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
      call rdone(irc,iopt,label,icomp,work(ita),toper)
      if (irc /= 0) call Error()
      call square(Work(ita),Work(isa),nb,1,nb)
      call dgemm_('N','N',nb,nb,nb,One,Work(isd),nb,Work(isa),nb,Zero,Work(isr),nb)
      do kdir=1,nb
        trace = trace+Work(isr+nb*(kdir-1)+kdir-1)
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
    write(u6,'(3E20.10)') (-amat(idir,jdir),jdir=1,3)
  end do
  write(u6,'(A)') '   ---------------------------------------------------------'
end do
call Add_Info('AMAT',AMAT,9,5)

call free_work(itd)
call free_work(isd)
call free_work(ita)
call free_work(isa)
call free_work(isr)
call clsone(irc,iopt)

return

contains

subroutine Error()
  write(u6,*) ' *** Error in subroutine cmp_hfc ***'
  write(u6,'(A,A)') '     Label = ',Label
  call Abend()
end subroutine Error

end subroutine cmp_hfc
