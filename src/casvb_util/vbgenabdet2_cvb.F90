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

subroutine vbgenabdet2_cvb(idetavb,idetbvb,iconfs,nconf,nconfion,ndetvb,nel,noe,nalf,nbet,norb,xalf,xbet,mingrph,maxgrph,inewocc, &
                           iaccm)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ndetvb, idetavb(ndetvb), idetbvb(ndetvb), nconf, noe, iconfs(noe,nconf), nel, nconfion(0:nel), nalf, nbet, &
                     norb, xalf(0:norb,0:nalf), xbet(0:norb,0:nbet), mingrph(0:norb), maxgrph(0:norb), inewocc(norb), iaccm(norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iaind, iaocc, iastr, ibind, ibocc, ibstr, iconf, ii, incr, incrdet, indx, ioff_nconf, ion, iorb, nalfsing, &
                     nbetsing, nelsing, nstring
logical(kind=iwp), parameter :: debug = .false.
integer(kind=iwp), external :: indget_cvb, mstacki_cvb

! Set xalf and xbet for indget:
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf,0)
  maxgrph(iorb) = min(iorb,nalf)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet,0)
  maxgrph(iorb) = min(iorb,nbet)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

! Transformation matrix VB structures  <->  determinants
incrdet = 0
ioff_nconf = 0
do ion=0,nel/2
  nelsing = nel-2*ion
  nalfsing = nalf-ion
  nbetsing = nbet-ion
  ! Skip if configurations are incompatible with Ms:
  if ((nalfsing >= 0) .and. (nbetsing >= 0) .and. (nelsing >= 0)) then

    ! Generate alpha/beta strings for singly occupied electrons:
    call icomb_cvb(nelsing,nalfsing,nstring)
    iastr = mstacki_cvb(nalfsing*nstring)
    ibstr = mstacki_cvb(nbetsing*nstring)
    call stringen_cvb(nelsing,nalfsing,iwork(iastr),iwork(ibstr))
    if (debug) then
      write(u6,*) ' ionicity=',ion,' nconf=',nconfion(ion)
      write(u6,*) ' check alpha strings :'
      do i=1,nstring
        write(u6,*) i,' => ',(iwork(ii+iastr-1+(i-1)*nalfsing),ii=1,nalfsing)
      end do
      write(u6,*) ' check beta strings :'
      do i=1,nstring
        write(u6,*) i,' => ',(iwork(ii+ibstr-1+(i-1)*nbetsing),ii=1,nbetsing)
      end do
    end if

    do iconf=ioff_nconf+1,ioff_nconf+nconfion(ion)
      call imove_cvb(iconfs(1,iconf),inewocc,norb)
      incr = 0
      do iorb=1,norb
        if (inewocc(iorb) == 1) then
          incr = incr+1
          iaccm(incr) = iorb
        end if
        inewocc(iorb) = max(0,inewocc(iorb)-1)
      end do

      ! Spin string loop:
      do indx=1,nstring

        ! Alpha index in full string space ...
        do i=1,nalfsing
          iaocc = iaccm(iwork(i+(indx-1)*nalfsing+iastr-1))
          inewocc(iaocc) = inewocc(iaocc)+1
        end do
        iaind = indget_cvb(inewocc,nalf,norb,xalf)
        do i=1,nalfsing
          iaocc = iaccm(iwork(i+(indx-1)*nalfsing+iastr-1))
          inewocc(iaocc) = inewocc(iaocc)-1
        end do

        ! Beta index in full string space ...
        do i=1,nbetsing
          ibocc = iaccm(iwork(i+(indx-1)*nbetsing+ibstr-1))
          inewocc(ibocc) = inewocc(ibocc)+1
        end do
        ibind = indget_cvb(inewocc,nbet,norb,xbet)
        do i=1,nbetsing
          ibocc = iaccm(iwork(i+(indx-1)*nbetsing+ibstr-1))
          inewocc(ibocc) = inewocc(ibocc)-1
        end do

        incrdet = incrdet+1
        idetavb(incrdet) = iaind
        idetbvb(incrdet) = ibind
      end do
    end do
    call mfreei_cvb(iastr)
  end if
  ioff_nconf = ioff_nconf+nconfion(ion)
end do
if (debug) then
  write(u6,*) ' idetavb='
  write(u6,'(10i6)') idetavb
  write(u6,*) ' idetbvb='
  write(u6,'(10i6)') idetbvb
end if

return

end subroutine vbgenabdet2_cvb
