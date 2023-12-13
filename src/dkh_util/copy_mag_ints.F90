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

subroutine copy_mag_ints(natoms)
!**************************************************************
!
! An interface for handling magnetic integrals from Gen1Int
!
!**************************************************************

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natoms
integer(kind=iwp) :: IDUM(1), nmag, irc, iopt, iCmp, icomp, toper, iat, jtri, ncomp, lu_one
character(len=8) :: Label
real(kind=wp), allocatable :: scrt(:)

! Copy magnetic integrals from ONEREL to ONEINT

lu_one = 2
iopt = 0
irc = -1
call OpnOne(irc,iopt,'ONEREL',lu_one)
if (irc /= 0) call Error()

! Primitive integrals stored on ONEREL
Label = 'MAGXP  1'
iOpt = ibset(0,sOpSiz)
iComp = 1
tOper = 255
! Integral dimensions
call iRdOne(irc,iOpt,Label,iComp,idum,tOper)
if (irc /= 0) call Error()
nmag = idum(1)
! Some scratch space
call mma_allocate(scrt,nmag+4,label='scratch')
iOpt = 0
nComp = 9
do iat=1,nAtoms
  do jtri=1,2
    if (jtri == 1) then
      write(Label,'(A,I3)') 'MAGXP',iat
    else
      write(Label,'(A,I3)') 'MAGPX',iat
    end if
    do iComp=1,nComp
      ! Read the primitives from ONEREL
      iCmp = iComp
      call RdOne(iRC,iOpt,Label,iCmp,scrt,tOper)
      if (iRC /= 0) call Error()
      ! Close ONEREL
      call ClsOne(iRC,iOpt)
      ! Open ONEINT
      call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
      if (iRC /= 0) call Error()
      ! Write the primitives to ONEINT ?
      call WrOne(iRC,iOpt,Label,iComp,scrt,tOper)
      call ClsOne(iRC,iOpt)
      call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
      if (iRC /= 0) call Error()
    end do
  end do
end do
call mma_deallocate(scrt)
call ClsOne(iRC,iOpt)

return

contains

subroutine Error()
  use Definitions, only: u6
  write(u6,*) ' *** Error in subroutine Copy_Mag_ints ***'
  write(u6,'(A,A)') '     Label = ',Label
  call Abend()
end subroutine Error

end subroutine copy_mag_ints
