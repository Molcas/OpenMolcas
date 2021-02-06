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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine creates a list of basis function labels.                *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine MkType()

use Genano_globals, only: MxLqn, symlab
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iFrom, iTo, l, m
integer(kind=iwp), parameter :: MxLst = MxLqn*(MxLqn*(MxLqn+6)+11)/6+1
character(len=8) :: CrtLst(0:MxLst), SphLst(0:MxLst)

!----------------------------------------------------------------------*
! Get labels from utility routine                                      *
!----------------------------------------------------------------------*
call Make_Labels(CrtLst,SphLst,MxLst,MxLqn)
!write(u6,*) 'Cartesian label'
!ind = 0
!do l=0,MxLqn
!  write(u6,'(15(1x,a4))') (CrtLst(ind+m),m=0,l*(l+3)/2)
!  ind = ind+l*(l+3)/2+1
!end do
!write(u6,*) 'Spherical label'
!ind = 0
!do i=0,MxLqn
!  do l=i,0,-2
!    write(u6,'(15(1x,a4))') (SphLst(ind+m),m=0,2*l)
!    ind = ind+2*l+1
!  end do
!end do
!----------------------------------------------------------------------*
! Copy labels                                                          *
!----------------------------------------------------------------------*
iFrom = 0
iTo = 1
do l=0,MxLqn
  if (l < 2) then
    do m=0,2*l
      symlab(iTo+m) = CrtLst(iFrom+m)
    end do
  else
    do m=0,2*l
      symlab(iTo+m) = SphLst(iFrom+m)
    end do
  end if
  iFrom = iFrom+l*(l+3)/2+1
  iTo = iTo+2*l+1
end do
!write(u6,*) 'Local labels'
!ind=1
!do l=0,MxLqn
!   write(u6,'(15(1x,a4))') (symlab(ind+m),m=0,2*l)
!   ind = ind+2*l+1
!end do

return

end subroutine Mktype
