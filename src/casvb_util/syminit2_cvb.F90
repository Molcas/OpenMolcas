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

implicit real*8(a-h,o-z)
logical found
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension symelm(norb,norb,nsyme), iorbrel(ndimrel)
dimension irels(2,*), relorb(norb,norb,*)
dimension north(norb), corth(norb,*)
dimension io(4,norbrel), iorder(norb,norbrel), iorbs(norb)
dimension a(norb,norb), b(norb,norb)
dimension rr(norb), ri(norb), vr(norb,norb), vi(norb,norb), intger(norb)
dimension ifxorb(norb), ifxstr(nfxvb), idelstr(nzrvb)
dimension iorts(2,nort), irots(2,ndrot), izeta(nsyme)
dimension dum(1)
save thresh
data thresh/1.d-8/

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
  if (iorb == jorb) goto 11
  call izero(iorbs,norb)
  iorbs(iorb) = 1
  iorbs(jorb) = 1
  ishift2 = 0
  do ijrel2=1,nijrel
    iorb2 = abs(iorbrel(1+ishift2))
    jorb2 = abs(iorbrel(2+ishift2))
    if ((ijrel2 == ijrel) .or. (iorb2 == jorb2)) goto 31
    if (iorbs(iorb2) == 1) then
      if (iorbs(jorb2) == 1) then
        ncount = 0
        do i=1,norb
          if (iorbs(i) == 1) then
            ncount = ncount+1
            intger(ncount) = i
          end if
        end do
        write(6,'(a,/,20i4)') ' Too many orbital relations involving orbitals :',(intger(ii),ii=1,ncount)
        write(6,'(a)') ' Please reduce number of ORBREL cards.'
        call abend_cvb()
      else
        iorbs(jorb2) = 1
      end if
    else if (iorbs(jorb2) == 1) then
      iorbs(iorb2) = 1
    end if
31  ishift2 = ishift2+3+iorbrel(3+ishift2)
  end do
11 ishift = ishift+3+iorbrel(3+ishift)
end do

! Diagonal orbital relations:
nciorth = 0
call izero(north,norb)
do iorb=1,norb
  ! Orbital conditions on IORB:
  call span0_cvb(norb,norb)
  ishift = 0
  do i=1,norbrel
    ior = iorbrel(1+ishift)
    jor = iorbrel(2+ishift)
    nrel = iorbrel(3+ishift)
    if ((iorb == ior) .and. (iorb == jor)) then
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
        write(6,*) ' Error in diagonalisation, IFAIL :',ifail
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
600 continue
! Sort relations and define generating orbitals:
do i=1,norb
  iorb = iorbs(i)
  do ii=1,nijrel
    if ((iorder(1,ii) == iorb) .or. (iorder(2,ii) == iorb)) then
      ! Has IORB already been generated from KORB?:
      do j=1,i-1
        korb = iorbs(j)
        if ((iorder(1,ii) == korb) .or. (iorder(2,ii) == korb)) goto 701
      end do
      if (iorder(1,ii) == iorb) then
        iorder(1,ii) = iorder(2,ii)
        iorder(2,ii) = iorb
      end if
      jorb = iorder(1,ii)
      ! JORB will be generated from IORB
      if (north(jorb) /= 0) then
        write(6,'(2(a,i4),a)') ' Attempting to generate orbital',jorb,' from orbital',iorb,'  ---'
        write(6,'(a,i4,a)') ' the orbital conditions for orbital',jorb,' cannot be enforced.'
        write(6,'(a)') ' Please reduce number of ORBREL cards.'
        call abend_cvb()
      end if
      found = .false.
      do jj=1,nijrel
        if ((jj /= ii) .and. ((iorder(1,jj) == jorb) .or. (iorder(2,jj) == jorb))) then
          ! Is JORB involved in any other orbital relations?:
          do j=1,i-1
            korb = iorbs(j)
            if ((iorder(1,jj) == korb) .or. (iorder(2,jj) == korb)) goto 900
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
        end if
900     continue
      end do
      if (found) goto 600
    end if
701 continue
  end do
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
    if (iorder(il,ijrel) == 0) goto 1200
    iorb = iorder(il,ijrel)
    j = i
1300 j = j+1
    jl = norb+2-j
    if (j == 1) jl = 2
    if (j == norb) jl = 1
    if (iorder(jl,ijrel) == 0) goto 1300
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
1200 continue
  end do
end do

return

end subroutine syminit2_cvb
!********************
!** Symmetrization **
!********************
