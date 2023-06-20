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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine dampF(r,RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)

use LEVEL_COMMON, only: bDS, cDS
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: r, RHOAB
integer(kind=iwp), intent(in) :: NCMM, MMLR(NCMM), IVSR, IDSTT
real(kind=wp), intent(out) :: DM(NCMM)
integer(kind=iwp) :: FIRST = 1, m, MM
real(kind=wp) :: br, XP, YP, ZK
real(kind=wp), save :: bpm(20,-2:0), cpm(20,-2:0)

!write(u6,*) 'Made it inside of dampF! IVSR=',IVSR
if (NCMM > 4) write(u6,*) 'IDSTT=',IDSTT
if (FIRST == 1) then
  do m=1,20
    bpm(m,:) = bDS(-2:0)/real(m,kind=wp)
    cpm(m,:) = cDS(-2:0)/sqrt(real(m,kind=wp))
  end do
  FIRST = 0
end if
br = RHOAB*r
!write(u6,*) 'NCMM=',NCMM
do m=1,NCMM
  MM = MMLR(m)
  XP = exp(-(bpm(MM,IVSR)+cpm(MM,IVSR)*br)*br)
  YP = One-XP
  ZK = MM-One
  DM(m) = YP**(MM-1)
  !... Actually ...  DM(m)= YP**(MM + IVSR/2)  :  set it up this way to
  !   avoid taking exponential of a logarithm for fractional powers (slow)
  select case (IVSR)
    case (-4)
      ZK = ZK-One
      DM(m) = DM(m)/YP
    case (-3)
      ZK = ZK-Half
      DM(m) = DM(m)/sqrt(YP)
    case (-1)
      ZK = ZK+Half
      DM(m) = DM(m)*sqrt(YP)
    case (0)
      ZK = MM
      DM(m) = DM(m)*YP
    case (-9)
  end select
end do

return

!600 format(/,' *** ERROR ***  For  IDSTT=',i3,'   IVSR=',i3,' no damping function is defined')
!602 format( /,' ***ERROR ***  RHOAB=', F7.4,'  yields an invalid Damping Function definition')

end subroutine dampF
