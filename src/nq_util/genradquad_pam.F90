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

subroutine GenRadQuad_PAM(iNQ,nR_Eff,mr,Alpha,Process,QuadR,nQuadR)

use nq_Info

implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "real.fh"
#include "debug.fh"
real*8 Alpha(2), QuadR(2,nQuadR)
real*8 mr(2), ln_rn
logical Process

! Last point at infinity is eliminated

if (Debug) write(6,*) 'New Algorithm (Malmqvist)'

! Reading of the data
Alpha_Min = Alpha(1)
Alpha_Max = Alpha(2)
l_Max = 2*int(mr(1))
if (Debug) write(6,*) 'l_Max=',l_Max
Relative_Max_Error = mr(2)
if (Debug) write(6,*) 'Relative_Max_Error=',Relative_Max_Error

! Compute an approximative R_D_0

Dr = Zero
h = Zero
do k=0,l_Max,l_Max-1
  R_D_0 = Relative_Max_Error/(10.0d0**k)
  Dr = -log10(R_D_0)

  ! Starting value of h

  h = One/(0.47d0*Dr+0.93d0)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute a correct h for the approximative R_D_0

  C1 = Four*sqrt(Two)*Pi
  C2 = Pi**2/Two
99 continue
  h_ = C2/(-log(Ten**(-Dr)*h/C1))
  if (abs(h_-h) > 1.0D-4) then
    h = h_
    Go To 99
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now find h from the correct R_D, i.e. from the highest
  ! angular momentum available.

  Dr = -log10(Relative_Max_Error)
98 continue
  h_ = C2/(-log(Ten**(-Dr)*(h/C1)*(h/Pi)**(dble(k)/Two)*(G((dble(k)+Three)/Two)/G(Three/Two))))

  if (Debug) write(6,*) 'h h_ ',h,h_
  if (abs(h_-h) > 1.0D-5) then
    h = h_
    Go To 98
  end if

  if (k == 0) h0 = h
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute table of R_Max as a function of l

do i=l_Max,0,-2
  D_m = -4.0d0
  if (l_Max == 4) D_m = -2.3d0
  if (l_Max == 2) D_m = -1.0d0
  if (l_Max == 0) D_m = 1.9d0
  if (l_Max == -2) D_m = 9.1d0

  ggg = (Two/(dble(i)+Three))*(D_m-log(One/Ten**(-Dr)))
  R_Max(i) = sqrt(exp(ggg)/Alpha_Max)
  if (Debug) then
    write(6,*) 'i        =',i
    write(6,*) 'l_Max    =',l_Max
    write(6,*) 'ggg      =',ggg
    write(6,*) 'R_Max(i) =',R_Max(i)
  end if
end do
if (Debug) write(6,*) 'h0,h=',h0,h

! For hybrid grid use R_Max for l=0 and h for l=l_max
! r1=R_Max(l_Max)
r1 = R_Max(0)
ln_rn = 1.7d0-log(Alpha_Min)/Two
rn = exp(ln_rn)
gamma = r1/(exp(h)-One)
n_High = int(log(rn/gamma+One)/h+One)
if (Debug) then
  write(6,*)
  write(6,*) 'r1,Alpha_Min    =',r1,Alpha_Min
  write(6,*) 'rn,Alpha_Max    =',rn,Alpha_MAx
  write(6,*) 'h,Dr,n_High     =',h,Dr,n_High
end if

! Store the radius and the associated weights

if (Debug) write(6,*) 'n_High',n_High
iR = 0
do k=0,n_High
  a = dble(k)*h
  rk = Gamma*(exp(a)-One)
  ! Note that the point at r=0 is eliminated
  if (rk /= Zero) then
    iR = iR+1
    if (Process) then
      QuadR(1,ir) = rk
      Correction = One

      ! Gregorious correction for points close to the nuclei

      if (k == 0) Correction = 46.d0/120.d0
      if (k == 1) Correction = 137.d0/120.d0
      if (k == 2) Correction = 118.d0/120.d0
      if (k == 3) Correction = 119.d0/120.d0

      QuadR(2,ir) = h*(rk+gamma)*Correction
      QuadR(2,ir) = QuadR(1,ir)**2*QuadR(2,ir)
    end if
  end if

end do
! Store the value of the maximum radius to which we should integrate
! for the partitionning and the number of effective radii.
nR_Eff = iR

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iNQ)

end subroutine GenRadQuad_PAM
