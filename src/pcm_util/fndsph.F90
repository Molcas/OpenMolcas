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

subroutine FndSph(NAt,ICharg,ToAng,C,IAt,ITypRad,NSphInp,Alpha,XSph,YSph,ZSph,Rad,NOrd,iPrint)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NAt, ICharg, IAt(*), ITypRad, NSphInp, iPrint
real(kind=wp), intent(in) :: ToAng, C(3,*)
real(kind=wp), intent(inout) :: Alpha
real(kind=wp), intent(_OUT_) :: XSph(*), YSph(*), ZSph(*), Rad(*)
integer(kind=iwp), intent(_OUT_) :: NOrd(*)
integer(kind=iwp) :: I
real(kind=wp), allocatable :: Chg(:)
real(kind=wp), external :: Pauling
#include "rctfld.fh"

! Assign GEPOL sphere positions and radii according to solute atoms nature
if (ITypRad == 1) then
  ! United Atom Topological Model (UATM) radii:
  call mma_allocate(Chg,NAt,label='Chg')
  Chg(:) = Zero
  call UATM(u6,ICharg,NAt,NSinit,ToAng,Rad,Alpha,C,IAt,NOrd,Chg,iPrint)
  call mma_deallocate(Chg)
elseif (ITypRad == 2) then
  ! Pauling radii on each atom:
  do I=1,NAt
    NOrd(I) = I
    Rad(I) = Pauling(IAt(I))
    Alpha = 1.2_wp
    NSinit = NAt
  end do
elseif (ITypRad == 3) then
  ! Sphere radii given in the input
  do I=1,NSphInp
    NOrd(I) = NOrdInp(I)
    Rad(I) = RadInp(I)
    Alpha = 1.2_wp
    NSinit = NSphInp
  end do
else
  write(u6,'(a)') 'Unrecognized radii type !'
  call Abend()
end if
if (((ITypRad == 2) .or. (ITypRad == 3)) .and. (iPrint > 5)) call PrtCav(u6,ITypRad,NSinit,NOrd,Alpha,Rad)

do I=1,NSinit
  XSph(I) = C(1,NOrd(I))
  YSph(I) = C(2,NOrd(I))
  ZSph(I) = C(3,NOrd(I))
  Rad(I) = Rad(I)*Alpha
end do

return

end subroutine FndSph
!====
subroutine PrtCav(IOut,ITyp,NS,NOrd,Alpha,Rad)
! Print out sphere radii for Pauling or input cavitites

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut, ITyp, NS, NOrd(*)
real(kind=wp) :: Alpha, Rad(*)
integer(kind=iwp) :: IS

write(iOut,*)
write(iOut,*)
write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
write(iOut,'(6X,A)') '================================'
if (ITyp == 2) write(iOut,'(6X,A)') 'Pauling radii'
if (ITyp == 3) write(iOut,'(6X,A)') 'Sphere radii from input'
write(iOut,*)
write(IOut,'(6X,A)') ' NOrd  Alpha  Radius'
do IS=1,NS
  write(IOut,'(6X,1X,I3,3X,F4.2,3X,F5.3)') NOrd(IS),Alpha,Rad(IS)
end do
write(IOut,'(6X,1X,78("-"))')
write(IOut,*)

return

end subroutine PrtCav
