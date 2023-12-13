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
! Copyright (C) 1990,1991, Roland Lindh                                *
!***********************************************************************

subroutine read_rysrw()
!***********************************************************************
!                                                                      *
! Object: to setup the coefficients for the Rys roots and weights.     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             September '90                                            *
!             Modified to DaFile February '91                          *
!***********************************************************************

use vRys_RW, only: Cff, ddx, iCffR, iCffW, iMap, ix0, Map, nMap, nMxRys, nx0, TMax, x0
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: i, io, iOff, iRys, lu_rysrw, mRys, nCff, nMap_Tot, nMem, nMem_Tot, nOrder, nx0_Tot
real(kind=wp) :: acc(size(iMap))
logical(kind=iwp) :: found_rysrw
character(len=*), parameter :: RYSRW_NAME = 'RYSRW'
integer(kind=iwp), external :: isFreeUnit

! Open file for data base

call f_Inquire(RYSRW_NAME,found_rysrw)
if (.not. found_rysrw) then
  call warningmessage(2,' the rysrw file does not exist.')
  call abend()
end if
lu_rysrw = isFreeUnit(22)
call molcas_open(lu_rysrw,RYSRW_NAME)

! Read initial data

io = 1
do while (io /= 0)
  read(lu_rysrw,*,iostat=io) mRys,nOrder
end do
if (mRys > size(iMap)) then
  call WarningMessage(2,' Database requires new code! Database and code are at incompatible levels!')
  call Abend()
end if
nMxRys = mRys
nCff = 2*(nOrder+1)
read(lu_rysrw,*) (Acc(i),i=1,mRys)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Reading tables for roots and weights of Rys poly.'
write(u6,*) ' Highest order is:',mRys
write(u6,*) ' Order of approximating polynomial:',nOrder
write(u6,*) ' Relative accuracy of computed values:',(Acc(i),i=1,mRys)
write(u6,*)
#endif

! Read value of T at which asymptotic formulas will be used

call mma_allocate(TMax,mRys,label='TMax')
read(lu_rysrw,*) (TMax(i),i=1,mRys)
#ifdef _DEBUGPRINT_
call RecPrt(' Tmax',' ',Tmax,mRys,1)
#endif

! Read increment of tables

call mma_allocate(ddx,mRys,label='ddx')
read(lu_rysrw,*) (ddx(i),i=1,mRys)
#ifdef _DEBUGPRINT_
call RecPrt(' ddx ',' ',ddx,mRys,1)
#endif

! Read size of map array

read(lu_rysrw,*) (nMap(i),i=1,mRys)
#ifdef _DEBUGPRINT_
write(u6,*) ' nMap=',nMap
#endif

! Read number of subranges

read(lu_rysrw,*) (nx0(i),i=1,mRys)
#ifdef _DEBUGPRINT_
write(u6,*) ' nx0=',nx0
#endif

! Read map array and x0 array for each order of Rys polynomials

nMap_Tot = 0
nx0_Tot = 0
do iRys=1,mRys
  iMap(iRys) = nMap_Tot+1
  nMap_Tot = nMap_Tot+nMap(iRys)
  ix0(iRys) = nx0_Tot+1
  nx0_Tot = nx0_Tot+nx0(iRys)
end do
call mma_allocate(Map,nMap_Tot,label='Map')
call mma_allocate(x0,nx0_Tot,label='x0')
do iRys=1,mRys
  iOff = iMap(iRys)-1
  read(lu_rysrw,*) (Map(i),i=iOff+1,iOff+nMap(iRys))

  iOff = ix0(iRys)-1
  read(lu_rysrw,*) (x0(i),i=iOff+1,iOff+nx0(iRys))
end do

! Allocate memory for coefficients

nMem_Tot = 0
do iRys=1,mRys
  iCffR(0,iRys) = nMem_Tot+1
  nMem = nx0(iRys)*iRys
  nMem_Tot = nMem_Tot+nCff*nMem
end do
call mma_allocate(Cff,nMem_Tot,label='Cff')
do iRys=1,mRys

  ! Read coefficients from file

  nMem = nx0(iRys)*iRys
  iCffR(1,iRys) = iCffR(0,iRys)+nMem
  iCffR(2,iRys) = iCffR(1,iRys)+nMem
  iCffR(3,iRys) = iCffR(2,iRys)+nMem
  iCffR(4,iRys) = iCffR(3,iRys)+nMem
  iCffR(5,iRys) = iCffR(4,iRys)+nMem
  iCffR(6,iRys) = iCffR(5,iRys)+nMem

  ICffW(0,iRys) = iCffR(6,iRys)+nMem
  iCffW(1,iRys) = iCffW(0,iRys)+nMem
  iCffW(2,iRys) = iCffW(1,iRys)+nMem
  iCffW(3,iRys) = iCffW(2,iRys)+nMem
  iCffW(4,iRys) = iCffW(3,iRys)+nMem
  iCffW(5,iRys) = iCffW(4,iRys)+nMem
  iCffW(6,iRys) = iCffW(5,iRys)+nMem

  iOff = iCffR(0,iRys)-1
  read(lu_rysrw,*) (Cff(i),i=iOff+1,iOff+nMem*nCff)

end do

close(lu_rysrw)

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unusued_real_array(Acc)
#endif

end subroutine read_rysrw
