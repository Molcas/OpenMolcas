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
! From: D.R. Yarkony, J. Phys. Chem. A 105 (2001) 6277-6293
! and: J. Chem. Theory Comput. 12 (2016) 3636-3653
!***********************************************************************

subroutine CI_Summary(Lu)

use Slapaf_Info, only: Gx, Gx0, NAC, Energy
use Slapaf_Parameters, only: CallLast, iter

implicit none
integer Lu, n, i
real*8, dimension(:), allocatable :: g, h, tmp
real*8 gg, hh, gh, sg, sh, dgh, deltagh, beta_ang, norm_g, norm_h, st, srel, shead, peaked, bif, aux
real*8, external :: dDot_
character(Len=2) LabA
character(Len=40) Description
#include "stdalloc.fh"
#include "real.fh"

n = 3*size(Gx,2)
call mma_Allocate(g,n)
call mma_Allocate(h,n)

! Compute the orthogonal branching vectors
! Note that d(E1-E0)/dx is stored, but we want to use d((E1-E0)/2)/dx
! (and forces instead of gradients)

gg = dDot_(n,Gx0(1,1,iter),1,Gx0(1,1,iter),1)*Quart
hh = dDot_(n,NAC(1,1,iter),1,NAC(1,1,iter),1)
gh = -dDot_(n,Gx0(1,1,iter),1,NAC(1,1,iter),1) !Factor 2 included
beta_ang = atan2(gh,gg-hh)*Half
call dCopy_(n,Gx0(1,1,iter),1,g,1)
call dScal_(n,-Half*cos(beta_ang),g,1)
call dAxpY_(n,sin(beta_ang),NAC(1,1,iter),1,g,1)
call dCopy_(n,NAC(1,1,iter),1,h,1)
call dScal_(n,cos(beta_ang),h,1)
call dAxpY_(n,Half*sin(beta_ang),Gx0(1,1,iter),1,h,1)
gg = dDot_(n,g,1,g,1)
hh = dDot_(n,h,1,h,1)
norm_g = sqrt(gg)
norm_h = sqrt(hh)
if (norm_g > 1.0D-12) then
  call dScal_(n,One/norm_g,g,1)
else
  call dCopy_(n,[Zero],0,g,1)
end if
if (norm_h > 1.0D-12) then
  call dScal_(n,One/norm_h,h,1)
else
  call dCopy_(n,[Zero],0,h,1)
end if
! Ensure that the asymmetry will be positive
! this fixes which vector is x and which is y
if (hh > gg) then
  call SwapVe(g,h,n)
  aux = gg
  gg = hh
  hh = aux
  aux = norm_g
  norm_g = norm_h
  norm_h = aux
end if
sg = -dDot_(n,Gx(1,1,iter),1,g,1)
sh = -dDot_(n,Gx(1,1,iter),1,h,1)
! Ensure that the tilt heading will be in the first quadrant
! this fixes the signs of the x and y vectors
if (sg < Zero) then
  sg = abs(sg)
  call DScal_(n,-One,g,1)
end if
if (sh < Zero) then
  sh = abs(sh)
  call DScal_(n,-One,h,1)
end if
st = sqrt(sg**2+sh**2)
dgh = sqrt((gg+hh)/Two)
LabA = ''
deltagh = gg-hh
if ((gg+hh) > 1.0D-12) then
  deltagh = deltagh/(gg+hh)
else
  LabA = ' *'
end if
if (dgh > 1.0D-12) then
  srel = st/dgh
else
  srel = Zero
end if
shead = atan2(sh,sg)

! peaked/sloped, bifurcating/single-path parameters

peaked = srel**2/(One-deltagh**2)*(One-deltagh*cos(Two*shead))
bif = ((One+deltagh)*cos(shead)**2)**(One/Three)
bif = bif+((One-deltagh)*sin(shead)**2)**(One/Three)
bif = (srel/(Two*deltagh))**(Two/Three)*bif
Description = ''
if (peaked < 1) then
  Description = trim(Description)//'peaked (P<1)'
else if (peaked > 1) then
  Description = trim(Description)//'sloped (P>1)'
else
  Description = trim(Description)//'* (P=1)'
end if
if (bif < 1) then
  Description = trim(Description)//' bifurcating (B<1)'
else if (bif > 1) then
  Description = trim(Description)//' single-path (B>1)'
else
  Description = trim(Description)//' * (B=1)'
end if

! Disable Last_Energy to prevent further rotations

CallLast = .false.

write(Lu,*)
call CollapseOutput(1,'Conical Intersection Characterization')
write(Lu,'(3X,A)') '-------------------------------------'
write(Lu,*)
write(Lu,*) 'See: J. Chem. Theory Comput. 12 (2016) 3636-3653'
write(Lu,*)
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Gradient difference','',Gx0(1,1,iter),n,1)
call RecPrt('Coupling vector','',NAC(1,1,iter),n,1)
call RecPrt('Average gradient','',Gx(1,1,iter),n,1)
write(Lu,100) 'Beta angle:',beta_ang
write(Lu,*)
#endif
write(Lu,100) 'Pitch (delta_gh):',dgh,' Eh/a0'
write(Lu,100) 'Asymmetry (Delta_gh):',deltagh,trim(LabA)
#ifdef _DEBUGPRINT_
write(Lu,100) 'Total tilt (s):',st,' Eh/a0'
#endif
write(Lu,100) 'Relative tilt (sigma=s/delta_gh):',srel
write(Lu,100) 'Tilt heading (theta_s):',shead
write(Lu,*)
write(Lu,101) 'P:',peaked
write(Lu,101) 'B:',bif
write(Lu,101) 'Type: '//trim(Description)
write(Lu,*)
write(Lu,*) 'Local linear representation:'
call mma_Allocate(tmp,n)
do i=1,size(Gx,2)
  tmp(0*size(Gx,2)+i) = g((i-1)*3+1)
  tmp(1*size(Gx,2)+i) = g((i-1)*3+2)
  tmp(2*size(Gx,2)+i) = g((i-1)*3+3)
end do
call RecPrt('Local x','',tmp,size(Gx,2),3)
do i=1,size(Gx,2)
  tmp(0*size(Gx,2)+i) = h((i-1)*3+1)
  tmp(1*size(Gx,2)+i) = h((i-1)*3+2)
  tmp(2*size(Gx,2)+i) = h((i-1)*3+3)
end do
call RecPrt('Local y','',tmp,size(Gx,2),3)
write(Lu,*)
write(Lu,110) Energy(iter),sg,sh
write(Lu,120) Two*dgh,deltagh
call mma_Deallocate(tmp)
call CollapseOutput(0,'Conical Intersection Characterization')
100 format(5X,A,T40,ES12.5,A)
101 format(5X,A,T11,ES12.5)
110 format(5X,'Average energy: ',F15.8,' + ',F12.8,'*x + ',F12.8,'*y')
120 format(5X,'Energy difference: ',F12.8,'*sqrt(r^2 + ',F12.8,'*t)',/,10X,'r^2 = x^2 + y^2',/,10X,'t = x^2 - y^2')

call mma_Deallocate(g)
call mma_Deallocate(h)

return

end subroutine CI_Summary
