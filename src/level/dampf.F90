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

use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NCMM, MMLR(NCMM), IVSR, IDSTT
real(kind=wp) :: r, RHOAB, DM(NCMM)
integer(kind=iwp) :: FIRST = 1, IDFF, m, MM
real(kind=wp) :: br, XP, YP, ZK
real(kind=wp), save :: bpm(20,-2:0), cpm(20,-2:0)
! what are these numbers?
real*8, parameter :: bDS(-4:0) = [2.50_wp,2.90_wp,3.3_wp,3.69_wp,3.95_wp], cDS(-4:0) = [0.468_wp,0.446_wp,0.423_wp,0.40_wp,0.39_wp]

!write(u6,*) 'Made it inside of dampF! IVSR=',IVSR
if (NCMM > 4) then
  write(u6,*) 'IDSTT=',IDSTT
end if
if (FIRST == 1) then
  do m=1,20
    do IDFF=-2,0
      bpm(m,IDFF) = bDS(IDFF)/real(m,kind=wp)
      cpm(m,IDFF) = cDS(IDFF)/sqrt(real(m,kind=wp))
    end do
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
  if (IVSR == -4) then
    ZK = ZK-One
    DM(m) = DM(m)/YP
  end if
  if (IVSR == -3) then
    ZK = ZK-Half
    DM(m) = DM(m)/sqrt(YP)
  end if
  if (IVSR == -1) then
    ZK = ZK+Half
    DM(m) = DM(m)*sqrt(YP)
  end if
  if (IVSR == 0) then
    ZK = MM
    DM(m) = DM(m)*YP
  end if
  if (IVSR == -9) then
  end if
end do

return

!600 format(/,' *** ERROR ***  For  IDSTT=',i3,'   IVSR=',i3,' no damping function is defined')
!602 format( /,' ***ERROR ***  RHOAB=', F7.4,'  yields an invalid Damping Function definition')

end subroutine dampF
