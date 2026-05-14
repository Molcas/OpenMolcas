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
!               1992, Per Ake Malmqvist                                *
!***********************************************************************

subroutine SetUpR(nRys)
!***********************************************************************
!                                                                      *
! Object: to setup the coefficients for the Rys roots and weights.     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             September '90                                            *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Cehmistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to DaFile February '91                          *
!                                                                      *
!     Added: Call to READAB, P.-A. Malmqvist March 1992, to set up     *
!             the tables needed to calculate large-order roots and     *
!             weights on request.                                      *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW, MaxHer
use vRys_RW, only: HerR2, HerW2, iHerR2, iHerW2
use abdata, only: read_abdata
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#ifdef _RYS_SCRATCH_
use RysScratch, only: SetAux
use Definitions, only: wp
#endif

implicit none
integer(kind=iwp), intent(in) :: nRys
integer(kind=iwp) :: iHer, iOffR, iRys, jRys, MemHer

if (allocated(iHerR2)) then
  call WarningMessage(2,'SetupR: Rys_Status is already active!')
  call Abend()
end if

#ifdef _RYS_SCRATCH_
call SetAux(1.0e-16_wp)
#endif

call Read_ABData()

call Read_RysRW()

! Set up the square of roots and the weights for Hermite polynomials
! We will only do this for the even numbered polynomials.

MemHer = nRys*(nRys+1)/2
call mma_allocate(iHerR2,nRys,label='iHerR2')
iHerR2(1) = 1
call mma_allocate(iHerW2,nRys,label='iHerW2')
iHerW2(1) = 1
call mma_allocate(HerR2,MemHer,label='HerR2')
call mma_allocate(HerW2,MemHer,label='HerW2')

if (2*nRys > MaxHer) then
  call WarningMessage(2,'SetupR: 2*nRys>MaxHer')
  call Abend()
end if
do iRys=1,nRys
  iHer = 2*iRys
  iOffR = (iRys*(iRys-1))/2
  iHerR2(iRys) = iHerR2(1)+iOffR
  iHerW2(iRys) = iHerW2(1)+iOffR
  do jRys=0,iRys-1
    HerR2(iHerR2(iRys)+jRys) = HerR(iHerR(iHer)+iRys+jRys)**2
    HerW2(iHerW2(iRys)+jRys) = HerW(iHerW(iHer)+iRys+jRys)
  end do
end do

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call TriPrt(' Hermite squared roots',' ',HerR2(iHerR2(1)),nRys)
call TriPrt(' Hermite weights      ',' ',HerW2(iHerW2(1)),nRys)
#endif

return

end subroutine SetUpR
