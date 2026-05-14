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

subroutine mksyminit_cvb()

use casvb_global, only: corth, idelstr, ifxorb, ifxstr, iorbrel, iorts, irels, irots, izeta, ndimrel, ndrot, nfxvb, nijrel, &
                        niorth, norb, norbrel, nort, north, nsyme, nzrvb, recinp, relorb, symelm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iaddr, icnt, ieig, ifail, ii, iior, ijrel, ijrel2, il, ioffs, iorb, iorb2, ir, irel, ishift, ishift2, j, &
                     jj, jl, jor, jorb, jorb2, korb, nciorth, ncount, nrel
real(kind=wp) :: dum(1)
logical(kind=iwp) :: found
integer(kind=iwp), allocatable :: intger(:), io(:,:), iorder(:,:), iorbs(:)
real(kind=wp), allocatable :: a(:,:), b(:,:), ri(:), rr(:), vi(:,:), vr(:,:)
real(kind=wp), parameter :: thresh = 1.0e-8_wp

! Restore arrays:
call rdioff_cvb(9,recinp,ioffs)
call rdis_cvb(iorbrel,ndimrel,recinp,ioffs)
call rdis_cvb(ifxorb,norb,recinp,ioffs)
call rdis_cvb(ifxstr,nfxvb,recinp,ioffs)
call rdis_cvb(idelstr,nzrvb,recinp,ioffs)
call rdis_cvb(iorts,2*nort,recinp,ioffs)
call rdis_cvb(irots,2*ndrot,recinp,ioffs)
call rdis_cvb(izeta,nsyme,recinp,ioffs)

call mma_allocate(iorbs,norb,label='iorbs')
call mma_allocate(intger,norb,label='intger')

! First check that minimum number of orbital relations has been given
! (no cycle should be complete):
ishift = 0
do ijrel=1,nijrel
  iorb = abs(iorbrel(1+ishift))
  jorb = abs(iorbrel(2+ishift))
  if (iorb /= jorb) then
    iorbs(:) = 0
    iorbs(iorb) = 1
    iorbs(jorb) = 1
    ishift2 = 0
    do ijrel2=1,nijrel
      iorb2 = abs(iorbrel(1+ishift2))
      jorb2 = abs(iorbrel(2+ishift2))
      if ((ijrel2 /= ijrel) .and. (iorb2 /= jorb2)) then
        if (iorbs(iorb2) == 1) then
          if (iorbs(jorb2) == 1) then
            ncount = 0
            do i=1,norb
              if (iorbs(i) == 1) then
                ncount = ncount+1
                intger(ncount) = i
              end if
            end do
            write(u6,'(a,/,20i4)') ' Too many orbital relations involving orbitals :',(intger(ii),ii=1,ncount)
            write(u6,'(a)') ' Please reduce number of ORBREL cards.'
            call abend_cvb()
          else
            iorbs(jorb2) = 1
          end if
        else if (iorbs(jorb2) == 1) then
          iorbs(iorb2) = 1
        end if
      end if
      ishift2 = ishift2+3+iorbrel(3+ishift2)
    end do
  end if
  ishift = ishift+3+iorbrel(3+ishift)
end do

call mma_allocate(a,norb,norb,label='a')
call mma_allocate(b,norb,norb,label='b')
call mma_allocate(rr,norb,label='rr')
call mma_allocate(ri,norb,label='ri')
call mma_allocate(vr,norb,norb,label='vr')
call mma_allocate(vi,norb,norb,label='vi')

! Diagonal orbital relations:
nciorth = 0
north(:) = 0
do iorb=1,norb
  ! Orbital conditions on IORB:
  call span0_cvb(norb,norb)
  ishift = 0
  do i=1,norbrel
    iior = iorbrel(1+ishift)
    jor = iorbrel(2+ishift)
    nrel = iorbrel(3+ishift)
    if ((iorb == iior) .and. (iorb == jor)) then
      call unitmat(b,norb)
      do ir=nrel,1,-1
        irel = iorbrel(ir+3+ishift)
        call mxatb_cvb(symelm(:,:,irel),b,norb,norb,norb,a)
        b(:,:) = a(:,:)
      end do
      ! Everything that hasn't got eigenvalue +1 will be orthogonalised away
      ! Unsymmetric diagonalisation:
      ifail = 0
      call f02agf(b,norb,norb,rr,ri,vr,norb,vi,norb,intger,ifail)
      if (ifail /= 0) then
        write(u6,*) ' Error in diagonalisation, IFAIL :',ifail
        call abend_cvb()
      end if
      do ieig=1,norb
        if ((abs(rr(ieig)-One) > thresh) .or. (abs(ri(ieig)) > thresh)) then
          vr(:,ieig) = vr(:,ieig)+vi(:,ieig)
          call span1_cvb(vr(1,ieig),1,dum,norb,0)
        end if
      end do
    end if
    ishift = ishift+3+nrel
  end do
  !if (plc_const) call rconstr_plc(iorb)
  call span2_cvb(corth(1,1+nciorth),north(iorb),dum,norb,0)
  nciorth = nciorth+north(iorb)
end do
niorth = nciorth

call mma_deallocate(rr)
call mma_deallocate(ri)
call mma_deallocate(vr)
call mma_deallocate(vi)
call mma_deallocate(intger)
call mma_allocate(io,4,norbrel,label='io')
call mma_allocate(iorder,norb,norbrel,label='iorder')

! Off-diagonal relations:
io(:,:) = 0
iorder(:,:) = 0
ijrel = 0
ishift = 0
do i=1,norbrel
  iorb = iorbrel(1+ishift)
  jorb = iorbrel(2+ishift)
  nrel = iorbrel(3+ishift)
  if (iorb /= jorb) then
    ijrel = ijrel+1
    io(1,ijrel) = iorb
    io(2,ijrel) = jorb
    io(3,ijrel) = nrel
    io(4,ijrel) = 3+ishift
    iorder(1,ijrel) = iorb
    iorder(2,ijrel) = jorb
    irel = iorbrel(4+ishift)
  end if
  ishift = ishift+3+nrel
end do
nijrel = ijrel

! Orbitals with constraints should be generating orbitals if possible:
icnt = 0
do iorb=1,norb
  if (north(iorb) > 0) then
    icnt = icnt+1
    iorbs(icnt) = iorb
  end if
end do
do iorb=1,norb
  if (north(iorb) == 0) then
    icnt = icnt+1
    iorbs(icnt) = iorb
  end if
end do
do
  ! Sort relations and define generating orbitals:
  loop1: do i=1,norb
    iorb = iorbs(i)
    loop2: do ii=1,nijrel
      if ((iorder(1,ii) /= iorb) .and. (iorder(2,ii) /= iorb)) cycle loop2
      ! Has IORB already been generated from KORB?:
      do j=1,i-1
        korb = iorbs(j)
        if ((iorder(1,ii) == korb) .or. (iorder(2,ii) == korb)) cycle loop2
      end do
      if (iorder(1,ii) == iorb) then
        iorder(1,ii) = iorder(2,ii)
        iorder(2,ii) = iorb
      end if
      jorb = iorder(1,ii)
      ! JORB will be generated from IORB
      if (north(jorb) /= 0) then
        write(u6,'(2(a,i4),a)') ' Attempting to generate orbital',jorb,' from orbital',iorb,'  ---'
        write(u6,'(a,i4,a)') ' the orbital conditions for orbital',jorb,' cannot be enforced.'
        write(u6,'(a)') ' Please reduce number of ORBREL cards.'
        call abend_cvb()
      end if
      found = .false.
      loop3: do jj=1,nijrel
        if ((jj == ii) .or. ((iorder(1,jj) /= jorb) .and. (iorder(2,jj) /= jorb))) cycle loop3
        ! Is JORB involved in any other orbital relations?:
        do j=1,i-1
          korb = iorbs(j)
          if ((iorder(1,jj) == korb) .or. (iorder(2,jj) == korb)) cycle loop3
        end do
        found = .true.
        if (iorder(1,jj) == jorb) then
          iorder(1,jj) = iorder(2,jj)
          iorder(2,jj) = jorb
        end if
        iorder(3,jj) = iorder(2,jj)
        iorder(4:,jj) = iorder(3:norb-1,ii)
        ! KORB will be generated from IORB (via JORB):
        iorder(2,jj) = iorder(2,ii)
      end do loop3
      if (found) exit loop1
    end do loop2
  end do loop1
  if (i > norb) exit
end do
! Generate transformation matrix for each relation:
irels(1:2,1:nijrel) = iorder(1:2,1:nijrel)
do ijrel=1,nijrel
  call unitmat(relorb(:,:,ijrel),norb)
  do i=1,norb-1
    il = norb+2-i
    if (i == 1) il = 2
    if (i == norb) il = 1
    if (iorder(il,ijrel) == 0) cycle
    iorb = iorder(il,ijrel)
    j = i
    do
      j = j+1
      jl = norb+2-j
      if (j == 1) jl = 2
      if (j == norb) jl = 1
      if (iorder(jl,ijrel) /= 0) exit
    end do
    jorb = iorder(jl,ijrel)
    do ii=1,nijrel
      if (((iorb == io(1,ii)) .or. (iorb == io(2,ii))) .and. ((jorb == io(1,ii)) .or. (jorb == io(2,ii)))) then
        nrel = io(3,ii)
        iaddr = io(4,ii)
        ! Operate right-to-left:
        do ir=nrel,1,-1
          irel = iorbrel(ir+iaddr)
          if (jorb == io(1,ii)) then
            call mxatb_cvb(symelm(:,:,irel),relorb(:,:,ijrel),norb,norb,norb,a)
          else
            call mxattb_cvb(symelm(:,:,irel),relorb(:,:,ijrel),norb,norb,norb,a)
          end if
          relorb(:,:,ijrel) = a(:,:)
        end do
      end if
    end do
  end do
end do

call mma_deallocate(io)
call mma_deallocate(iorder)
call mma_deallocate(iorbs)
call mma_deallocate(a)
call mma_deallocate(b)

return

end subroutine mksyminit_cvb
