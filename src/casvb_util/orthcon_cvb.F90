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

subroutine orthcon_cvb(ipairs,mxortl,mxpair)

use casvb_global, only: mxorb_cvb, nort
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mxpair, mxortl
integer(kind=iwp), intent(out) :: ipairs(2,mxpair)
integer(kind=iwp), parameter :: ncmp = 4, nstrin = 7, mxgroup = 40
integer(kind=iwp) :: i, igrp, io, ior1, ior2, ipar, isp, istr, j, jo, jor1, jor2, jsp, ngrp, npairs, nread, nsp
character(len=3) :: glabel(mxgroup)
integer(kind=iwp), allocatable :: igroups(:,:), iorthlst(:), ipair(:,:), ngroup(:)
character(len=*), parameter :: string(nstrin) = ['GROUP   ','ORTH    ','PAIRS   ','STRONG  ','FULL    ','END     ','ENDORTHC']

call mma_allocate(ipair,mxorb_cvb,mxorb_cvb,label='ipair')
call mma_allocate(igroups,mxorb_cvb,mxgroup,label='igroups')
call mma_allocate(ngroup,mxgroup,label='ngroup')

ipair(:,:) = 0
ngrp = 0
do
  call fstring_cvb(string,nstrin,istr,ncmp,2)
  if (istr == 1) then
    ! 'GROUP'
    ngrp = ngrp+1
    if (ngrp > mxgroup) then
      write(u6,*) ' Too many GROUP keywords in input!',mxgroup
      call abend_cvb()
    end if
    glabel(ngrp) = ' '
    call string_cvb(glabel(ngrp),1,nread,1)
    if ((glabel(ngrp)(1:1) < 'A') .or. (glabel(ngrp)(1:1) > 'Z')) then
      write(u6,*) ' Group label must begin with a character A-Z: ',glabel(ngrp)
      call abend_cvb()
    end if
    call int_cvb(igroups(1,ngrp),mxorb_cvb,ngroup(ngrp),0)
    if (ngroup(ngrp) == -1) then
      write(u6,*) ' Too many elements for group ',glabel(ngrp)
      call abend_cvb()
    end if
    do i=1,ngroup(ngrp)
      if ((igroups(i,ngrp) < 1) .or. (igroups(i,ngrp) > mxorb_cvb)) then
        write(u6,*) ' Illegal orbital number in group ',glabel(ngrp),' :',igroups(i,ngrp)
        call abend_cvb()
      end if
    end do
    do i=1,ngrp-1
      if (glabel(ngrp) == glabel(i)) then
        write(u6,*) ' Repeated label in GROUP keywords : ',glabel(ngrp)
        call abend_cvb()
      end if
    end do
  else if (istr == 2) then
    ! 'ORTH'
    nsp = 0
    call mma_allocate(iorthlst,mxortl,label='iorthlst')
    do
      call int_cvb(iorthlst(1+nsp),mxortl-nsp,nread,0)
      nsp = nsp+nread
      if (mxortl-nsp <= 0) exit
      call fstring_cvb(glabel,ngrp,igrp,3,0)
      if (igrp <= 0) exit
      nsp = nsp+1
      iorthlst(nsp) = -igrp
      if (mxortl-nsp <= 0) exit
    end do
    do isp=1,nsp
      do jsp=isp+1,nsp
        ior1 = iorthlst(isp)
        jor1 = iorthlst(jsp)
        if ((ior1 > 0) .and. (jor1 > 0)) then
          ipair(ior1,jor1) = 1
        else if ((ior1 > 0) .and. (jor1 < 0)) then
          jor2 = -jor1
          do jo=1,ngroup(jor2)
            ipair(ior1,igroups(jo,jor2)) = 1
          end do
        else if ((ior1 < 0) .and. (jor1 > 0)) then
          ior2 = -ior1
          do io=1,ngroup(ior2)
            ipair(jor1,igroups(io,ior2)) = 1
          end do
        else if ((ior1 < 0) .and. (jor1 < 0)) then
          ior2 = -ior1
          jor2 = -jor1
          do io=1,ngroup(ior2)
            do jo=1,ngroup(jor2)
              ipair(igroups(io,ior2),igroups(jo,jor2)) = 1
            end do
          end do
        end if
      end do
    end do
    call mma_deallocate(iorthlst)
  else if (istr == 3) then
    ! 'PAIRS'
    nsp = 0
    do
      call int_cvb(ipairs(1+nsp,1),2*mxpair-nsp,nread,0)
      nsp = nsp+nread
      if (2*mxpair-nsp <= 0) exit
      call fstring_cvb(glabel,ngrp,igrp,3,0)
      if (igrp <= 0) exit
      nsp = nsp+1
      ipairs(nsp,1) = -igrp
      if (mxortl-nsp <= 0) exit
    end do
    if (mod(nsp,2) == 1) then
      write(u6,*) ' Odd number of orthogonalization numbers in PAIRS!'
      call abend_cvb()
    end if
    npairs = nsp/2
    do ipar=1,npairs
      ior1 = ipairs(1,ipar)
      jor1 = ipairs(2,ipar)
      if ((ior1 > 0) .and. (jor1 > 0)) then
        ipair(ior1,jor1) = 1
      else if ((ior1 > 0) .and. (jor1 < 0)) then
        jor2 = -jor1
        do jo=1,ngroup(jor2)
          ipair(ior1,igroups(jo,jor2)) = 1
        end do
      else if ((ior1 < 0) .and. (jor1 > 0)) then
        ior2 = -ior1
        do io=1,ngroup(ior2)
          ipair(jor1,igroups(io,ior2)) = 1
        end do
      else if ((ior1 < 0) .and. (jor1 < 0)) then
        ior2 = -ior1
        jor2 = -jor1
        do io=1,ngroup(ior2)
          do jo=1,ngroup(jor2)
            ipair(igroups(io,ior2),igroups(jo,jor2)) = 1
          end do
        end do
      end if
    end do
  else if (istr == 4) then
    ! 'STRONG'
    do i=1,mxorb_cvb
      do j=i+1,mxorb_cvb
        if (.not. ((mod(i,2) == 1) .and. (j == i+1))) ipair(i,j) = 1
      end do
    end do
  else if (istr == 5) then
    ! 'FULL'
    do i=1,mxorb_cvb
      ipair(i,i+1:) = 1
    end do
  end if
  ! 'END' , 'ENDORTHC' or unrecognized keyword -- end of ORTHCON input:
  if ((istr == 6) .or. (istr == 7) .or. (istr == 0)) exit
end do
ipairs(:,:) = 0
nort = 0
do i=1,mxorb_cvb
  do j=i+1,mxorb_cvb
    if ((ipair(i,j) == 1) .or. (ipair(j,i) == 1)) then
      nort = nort+1
      ipairs(1,nort) = i
      ipairs(2,nort) = j
    end if
  end do
end do

call mma_deallocate(ipair)
call mma_deallocate(igroups)
call mma_deallocate(ngroup)

return

end subroutine orthcon_cvb
