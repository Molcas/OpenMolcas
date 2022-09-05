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

subroutine GenRadQuad_PAM(nR_Eff,mr,Alpha,Process,QuadR,nQuadR)

use nq_Info, only: R_Max
use Constants, only: Zero, One, Two, Three, Four, Ten, Pi
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: nR_Eff
real(kind=wp), intent(in) :: mr(2), Alpha(2)
integer(kind=iwp), intent(in) :: nQuadR
real(kind=wp), intent(out) :: QuadR(2,nQuadR)
logical(kind=iwp) :: Process
integer(kind=iwp) :: i, iR, k, l_Max, n_High
real(kind=wp) :: a, Alpha_Max, Alpha_Min, C1, C2, Correction, D_m, Dr, Gmma, ggg, h, h_, ln_rn, r1, R_D_0, Relative_Max_Error, rk, &
                 rn
real(kind=wp), external :: G

! Last point at infinity is eliminated

#ifdef _DEBUGPRINT_
write(u6,*) 'New Algorithm (Malmqvist)'
#endif

! Reading of the data
Alpha_Min = Alpha(1)
Alpha_Max = Alpha(2)
l_Max = 2*int(mr(1))
#ifdef _DEBUGPRINT_
write(u6,*) 'l_Max=',l_Max
#endif
Relative_Max_Error = mr(2)
#ifdef _DEBUGPRINT_
write(u6,*) 'Relative_Max_Error=',Relative_Max_Error
#endif

! Compute an approximative R_D_0

Dr = Zero
h = Zero
do k=0,l_Max,l_Max-1
  R_D_0 = Relative_Max_Error/(Ten**k)
  Dr = -log10(R_D_0)

  ! Starting value of h

  h = One/(0.47_wp*Dr+0.93_wp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute a correct h for the approximative R_D_0

  C1 = Four*sqrt(Two)*Pi
  C2 = Pi**2/Two
  do
    h_ = C2/(-log(Ten**(-Dr)*h/C1))
    if (abs(h_-h) <= 1.0e-4_wp) exit
    h = h_
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now find h from the correct R_D, i.e. from the highest
  ! angular momentum available.

  Dr = -log10(Relative_Max_Error)
  do
    h_ = C2/(-log(Ten**(-Dr)*(h/C1)*(h/Pi)**(real(k,kind=wp)/Two)*(G((real(k,kind=wp)+Three)/Two)/G(Three/Two))))

#   ifdef _DEBUGPRINT_
    write(u6,*) 'h h_ ',h,h_
#   endif
    if (abs(h_-h) <= 1.0e-5_wp) exit
    h = h_
  end do

# ifdef _DEBUGPRINT_
  if (k == 0) h0 = h
# endif
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute table of R_Max as a function of l

do i=l_Max,0,-2
  D_m = -Four
  if (l_Max == 4) D_m = -2.3_wp
  if (l_Max == 2) D_m = -One
  if (l_Max == 0) D_m = 1.9_wp
  if (l_Max == -2) D_m = 9.1_wp

  ggg = (Two/(real(i,kind=wp)+Three))*(D_m-log(One/Ten**(-Dr)))
  R_Max(i) = sqrt(exp(ggg)/Alpha_Max)
# ifdef _DEBUGPRINT_
  write(u6,*) 'i        =',i
  write(u6,*) 'l_Max    =',l_Max
  write(u6,*) 'ggg      =',ggg
  write(u6,*) 'R_Max(i) =',R_Max(i)
# endif
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'h0,h=',h0,h
#endif

! For hybrid grid use R_Max for l=0 and h for l=l_max
! r1=R_Max(l_Max)
r1 = R_Max(0)
ln_rn = 1.7_wp-log(Alpha_Min)/Two
rn = exp(ln_rn)
Gmma = r1/(exp(h)-One)
n_High = int(log(rn/Gmma+One)/h+One)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'r1,Alpha_Min    =',r1,Alpha_Min
write(u6,*) 'rn,Alpha_Max    =',rn,Alpha_MAx
write(u6,*) 'h,Dr,n_High     =',h,Dr,n_High
#endif

! Store the radius and the associated weights

iR = 0
do k=0,n_High
  a = real(k,kind=wp)*h
  rk = Gmma*(exp(a)-One)
  ! Note that the point at r=0 is eliminated
  if (rk /= Zero) then
    iR = iR+1
    if (Process) then
      QuadR(1,ir) = rk
      Correction = One

      ! Gregorious correction for points close to the nuclei

      if (k == 0) Correction = 46.0_wp/120.0_wp
      if (k == 1) Correction = 137.0_wp/120.0_wp
      if (k == 2) Correction = 118.0_wp/120.0_wp
      if (k == 3) Correction = 119.0_wp/120.0_wp

      QuadR(2,ir) = h*(rk+Gmma)*Correction
      QuadR(2,ir) = QuadR(1,ir)**2*QuadR(2,ir)
    end if
  end if

end do
! Store the value of the maximum radius to which we should integrate
! for the partitionning and the number of effective radii.
nR_Eff = iR

return

end subroutine GenRadQuad_PAM
