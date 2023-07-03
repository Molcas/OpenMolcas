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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine Dihedr(Lbls,xyz,mCentr,rtrnc,Max_Center)
!***********************************************************************
!                                                                      *
! Object: to compute dihedral angles from a list of coordinates.       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "Molcas.fh"
real*8 xyz(3,mCentr), Bt(3,4), Coor(3,4)
character*(LenIn) Lbls(mCentr)
character*8 Label
logical type
dimension Dummy(1)

Lu = 6
Label = ' '
if (mCentr > Max_Center) Go To 99

Thr = 1.0D-12
type = .false.
do ic=1,mCentr
  x2 = xyz(1,ic)
  y2 = xyz(2,ic)
  z2 = xyz(3,ic)
  call dcopy_(3,xyz(1,ic),1,Coor(1,2),1)
  do jc=1,mCentr
    if (jc == ic) Go To 453
    x3 = xyz(1,jc)
    y3 = xyz(2,jc)
    z3 = xyz(3,jc)
    r2 = sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
    if ((r2 > rtrnc) .or. (r2 == Zero)) Go To 453
    !write(Lu,*)
    call dcopy_(3,xyz(1,jc),1,Coor(1,3),1)
    do kc=1,mCentr
      if (kc == ic) Go To 454
      if (kc == jc) Go To 454
      x1 = xyz(1,kc)
      y1 = xyz(2,kc)
      z1 = xyz(3,kc)
      r1 = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      if ((r1 > rtrnc) .or. (r1 == Zero)) Go To 454
      arg = ((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+(z1-z2)*(z3-z2))/(r1*r2)
      if (abs(arg) > One) arg = sign(One,arg)
      if (One-abs(arg) < Thr) Go To 454
      Phi1 = 180.d0*acos(arg)/Pi
      x12 = (y2-y1)*(z3-z2)-(y3-y2)*(z2-z1)
      y12 = (z2-z1)*(x3-x2)-(z3-z2)*(x2-x1)
      z12 = (x2-x1)*(y3-y2)-(x3-x2)*(y2-y1)
      r12 = sqrt(x12**2+y12**2+z12**2)
      if (r12 == Zero) Go To 454
      call dcopy_(3,xyz(1,kc),1,Coor(1,1),1)
      do lc=kc+1,mCentr
        if (lc == ic) Go To 455
        if (lc == jc) Go To 455
        if (lc == kc) Go To 455
        x4 = xyz(1,lc)
        y4 = xyz(2,lc)
        z4 = xyz(3,lc)
        r3 = sqrt((x4-x3)**2+(y4-y3)**2+(z4-z3)**2)
        if ((r3 > rtrnc) .or. (r3 == Zero)) Go To 455
        arg = ((x2-x3)*(x4-x3)+(y2-y3)*(y4-y3)+(z2-z3)*(z4-z3))/(r2*r3)
        if (abs(arg) > One) arg = sign(One,arg)
        if (One-abs(arg) < Thr) Go To 455
        Phi2 = 180.d0*acos(arg)/Pi
        x23 = (y3-y2)*(z4-z3)-(y4-y3)*(z3-z2)
        y23 = (z3-z2)*(x4-x3)-(z4-z3)*(x3-x2)
        z23 = (x3-x2)*(y4-y3)-(x4-x3)*(y3-y2)
        r23 = sqrt(x23**2+y23**2+z23**2)
        if (r23 == Zero) Go To 455
        call dcopy_(3,xyz(1,lc),1,Coor(1,4),1)
        !arg = (x12*x23+y12*y23+z12*z23)/(r12*r23)
        !if (abs(arg) > One) arg = sign(One,arg)
        !Phi12 = 180.D0*acos(arg)/Pi
        call Trsn(Coor,4,Tau,Bt,.false.,.false.,Label,Dummy,.false.)
        Phi12 = 180.0D+00*Tau/Pi
        if (.not. type) then
          type = .true.
          write(Lu,*)
          write(Lu,'(10X,A)') ' ***************************************************************'
          write(Lu,'(10X,A)') ' *              Valence Dihedral Angles / Degree               *'
          write(Lu,'(10X,A)') ' ***************************************************************'
          write(Lu,'(7X,A)') '             Atom centers                       Phi1     Phi2     Theta '
        end if
        write(Lu,'(10X,4(I2,1X,A,2X),1X,3(F7.2,2X))')kc,Lbls(kc),ic,Lbls(ic),jc,Lbls(jc),lc,Lbls(lc),Phi1,Phi2,Phi12
455     continue
      end do
454   continue
    end do
453 continue
  end do
end do

99 continue

return

end subroutine Dihedr
