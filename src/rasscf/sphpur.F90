!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SPHPUR(CMO)

use define_af, only: AngTp, iTabMx
use rasscf_global, only: BName, IXSYM
use general_data, only: NBAS, NORB, NSYM
use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp) :: I, IB, IBAS, IBASES, ICMOES, IO, IORB, IORBES, ISSLAB, ISYM, ITP, L, LCOUNT, LEXIST, LMX, MNL, MXL, NB, &
                     NBTOT, NO, NONZ
real(kind=wp) :: WGT, WGTLQN(0:9), WMX
logical(kind=iwp) :: IFTEST
character :: LCHAR
integer(kind=iwp), allocatable :: LQN(:)

! Set IFTEST=.true. to get supsym input generated in the output
! for further use, or for testing.
IFTEST = .false.
if (IFTEST) write(u6,*) 'SUPSYM'

! Set up array with angular quant num for each basis function:
NBTOT = sum(NBAS(1:NSYM))
call mma_allocate(LQN,NBTOT,Label='LQN')
do IBAS=1,NBTOT
  LCHAR = BName(IBAS)(LenIn+3:LenIn+3)
  L = -999999
  do ITP=0,ITABMX
    if (LCHAR == ANGTP(ITP)) L = ITP
  end do
  LQN(IBAS) = L
end do
ICMOES = 0
IBASES = 0
IORBES = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  NO = NORB(ISYM)
  if (NO /= 0) then
    do IO=1,NO
      IORB = IORBES+IO
      WGTLQN(:) = Zero
      do IB=1,NB
        IBAS = IBASES+IB
        L = LQN(IBAS)
        WGT = CMO(ICMOES+IB+NB*(IO-1))**2
        WGTLQN(L) = WGTLQN(L)+WGT
      end do
      LMX = 0
      WMX = WGTLQN(0)
      do L=0,9
        if (WGTLQN(L) > WMX) then
          LMX = L
          WMX = WGTLQN(L)
        end if
      end do
      IXSYM(IORB) = LMX
    end do
    ! We have now a provisional IXSYM array. How many different
    ! L values appear in it?
    MNL = 9
    MXL = 0
    NONZ = 0
    do L=0,9
      LEXIST = 0
      do IO=1,NO
        IORB = IORBES+IO
        if (L == IXSYM(IORB)) then
          LEXIST = 1
          MNL = min(L,MNL)
          MXL = max(L,MXL)
          exit
        end if
      end do
      NONZ = NONZ+LEXIST
    end do
    ! There are NONZ different values, so we want NONZ-1 special supsym
    ! labels for this symmetry. Reuse IWORK(LLQN) for orbital numbers:
    ! This will be the supsym label:
    if (IFTEST) write(u6,*) NONZ-1
    ISSLAB = 0
    do L=MNL,MXL
      LCOUNT = 0
      do IO=1,NO
        IORB = IORBES+IO
        if (L == IXSYM(IORB)) then
          LCOUNT = LCOUNT+1
          LQN(LCOUNT) = IO
        end if
      end do
      if (LCOUNT > 0) then
        ! Replace provisional IXSYM value with correct label:
        do IO=1,NO
          IORB = IORBES+IO
          if (IXSYM(IORB) == L) IXSYM(IORB) = ISSLAB
        end do
        ! Lowest L = label zero = do not specify in input:
        if (IFTEST .and. (ISSLAB > 0)) write(u6,'(1x,I3,16I5,(/,5X,16I5))') LCOUNT,(LQN(i),i=1,LCOUNT)
        ISSLAB = ISSLAB+1
      end if
    end do

    ICMOES = ICMOES+NO*NB
    IORBES = IORBES+NO
  end if
  IBASES = IBASES+NB
end do
call mma_deallocate(LQN)

end subroutine SPHPUR
