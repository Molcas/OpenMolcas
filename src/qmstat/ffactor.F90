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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************

!----------------------------------------------------------------------*
! A subroutine that computes those darn f-factors. They are defined    *
! in equation (2.4) in doi:10.1143/JPSJ.21.2313. As can be seen from   *
! that equation, the computation of the f-factors is actually a matter *
! of using the binomial theorem. This is what we do below and to make  *
! the computation efficient the expression (2.4) is written as a       *
! succint double sum.                                                  *
!----------------------------------------------------------------------*
subroutine fFactor(loneX,ltwoX,lsumX,loneY,ltwoY,lsumY,loneZ,ltwoZ,lsumZ,PAxyz,PBxyz,FactorX,FactorY,FactorZ)

use qmstat_global, only: MxAngqNr
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: loneX, ltwoX, lsumX, loneY, ltwoY, lsumY, loneZ, ltwoZ, lsumZ
real(kind=wp), intent(in) :: PAxyz(3), PBxyz(3)
real(kind=wp), intent(out) :: FactorX(2*MxAngqNr+1), FactorY(2*MxAngqNr+1), FactorZ(2*MxAngqNr+1)
integer(kind=iwp) :: i, ia, iLowB, iUpB
real(kind=wp) :: fff1, fff2, PAraise, PBraise
integer(kind=iwp), external :: NoverP_Q

do ia=0,lsumX !We use unrolled loops with regard to x,y and z therefore, here we start with the x-factors.
  fff2 = 0
  ! These lower and upper bounds have to do
  ! with the allowed numbers in the binomial coefficients.
  iLowB = max(0,ia-ltwoX)
  iUpB = min(ia,loneX)
  do i=iLowB,iUpB
    fff1 = NoverP_Q(loneX,i)*NoverP_Q(ltwoX,ia-i)
    if (i /= 0) then  !This is needed for some compilers (NAG_64)
      PAraise = PAxyz(1)**i
    else
      PAraise = One
    end if
    if (ia-i /= 0) then
      PBraise = PBxyz(1)**(ia-i)
    else
      PBraise = One
    end if
    fff2 = fff2+fff1*PAraise*PBraise
  end do
  FactorX(lsumX-ia+1) = fff2
end do
do ia=0,lsumY !y-factors.
  fff2 = 0
  iLowB = max(0,ia-ltwoY)
  iUpB = min(ia,loneY)
  do i=iLowB,iUpB
    fff1 = NoverP_Q(loneY,i)*NoverP_Q(ltwoY,ia-i)
    if (i /= 0) then !This is needed for some compilers (NAG_64)
      PAraise = PAxyz(2)**i
    else
      PAraise = One
    end if
    if (ia-i /= 0) then
      PBraise = PBxyz(2)**(ia-i)
    else
      PBraise = One
    end if
    fff2 = fff2+fff1*PAraise*PBraise
  end do
  FactorY(lsumY-ia+1) = fff2
end do
do ia=0,lsumZ !z-factorz.
  fff2 = 0
  iLowB = max(0,ia-ltwoZ)
  iUpB = min(ia,loneZ)
  do i=iLowB,iUpB
    fff1 = NoverP_Q(loneZ,i)*NoverP_Q(ltwoZ,ia-i)
    if (i /= 0) then !This is needed for some compilers (NAG_64)
      PAraise = PAxyz(3)**i
    else
      PAraise = One
    end if
    if (ia-i /= 0) then
      PBraise = PBxyz(3)**(ia-i)
    else
      PBraise = One
    end if
    fff2 = fff2+fff1*PAraise*PBraise
  end do
  Factorz(lsumZ-ia+1) = fff2
end do

return

end subroutine fFactor
