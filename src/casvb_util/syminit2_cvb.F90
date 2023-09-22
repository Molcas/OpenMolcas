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

subroutine syminit2_cvb(symelm,iorbrel,north,corth,irels,relorb,io,iorder,iorbs,a,b,rr,ri,vr,vi,intger,ifxorb,ifxstr,idelstr, &
                        iorts,irots,izeta)

use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: symelm(norb,norb,nsyme), corth(norb,*), relorb(norb,norb,*), a(norb,norb), b(norb,norb), rr(norb), ri(norb), &
                 vr(norb,norb), vi(norb,norb)
integer(kind=iwp) :: iorbrel(ndimrel), north(norb), irels(2,*), io(4,norbrel), iorder(norb,norbrel), iorbs(norb), intger(norb), &
                     ifxorb(norb), ifxstr(nfxvb), idelstr(nzrvb), iorts(2,nort), irots(2,ndrot), izeta(nsyme)
#include "files_cvb.fh"
integer(kind=iwp) :: i, iaddr, icnt, ieig, ifail, ii, iior, ijrel, ijrel2, il, ioffs, iorb, iorb2, ir, irel, ishift, ishift2, j, &
                     jj, jl, jor, jorb, jorb2, korb, nciorth, ncount, nrel
real(kind=wp) :: dum(1)
logical(kind=iwp) :: found
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

! First check that minimum number of orbital relations has been given
! (no cycle should be complete):
ishift = 0
do ijrel=1,nijrel
  iorb = abs(iorbrel(1+ishift))
  jorb = abs(iorbrel(2+ishift))
  if (iorb /= jorb) then
    call izero(iorbs,norb)
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

! Diagonal orbital relations:
nciorth = 0
call izero(north,norb)
do iorb=1,norb
  ! Orbital conditions on IORB:
  call span0_cvb(norb,norb)
  ishift = 0
  do i=1,norbrel
    iior = iorbrel(1+ishift)
    jor = iorbrel(2+ishift)
    nrel = iorbrel(3+ishift)
    if ((iorb == iior) .and. (iorb == jor)) then
      call mxunit_cvb(b,norb)
      do ir=nrel,1,-1
        irel = iorbrel(ir+3+ishift)
        call mxatb_cvb(symelm(1,1,irel),b,norb,norb,norb,a)
        call fmove_cvb(a,b,norb*norb)
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
        if ((abs(rr(ieig)-one) > thresh) .or. (abs(ri(ieig)) > thresh)) then
          call addvec(vr(1,ieig),vr(1,ieig),vi(1,ieig),norb)
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

! Off-diagonal relations:
call izero(io,2*norbrel)
call izero(iorder,norb*norbrel)
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
        call imove_cvb(iorder(3,ii),iorder(4,jj),norb-3)
        ! KORB will be generated from IORB (via JORB):
        iorder(2,jj) = iorder(2,ii)
      end do loop3
      if (found) exit loop1
    end do loop2
  end do loop1
  if (i > norb) exit
end do
! Generate transformation matrix for each relation:
do i=1,nijrel
  irels(1,i) = iorder(1,i)
  irels(2,i) = iorder(2,i)
end do
do ijrel=1,nijrel
  call mxunit_cvb(relorb(1,1,ijrel),norb)
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
            call mxatb_cvb(symelm(1,1,irel),relorb(1,1,ijrel),norb,norb,norb,a)
          else
            call mxattb_cvb(symelm(1,1,irel),relorb(1,1,ijrel),norb,norb,norb,a)
          end if
          call fmove_cvb(a,relorb(1,1,ijrel),norb*norb)
        end do
      end if
    end do
  end do
end do

return

end subroutine syminit2_cvb
