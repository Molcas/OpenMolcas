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

subroutine FndSph(NAt,ICharg,C,IAt,ITypRad,NSphInp,Alpha,XSph,YSph,ZSph,Rad,NOrd,m,iPrint)

use Solvent_Data, only: Pauling
use rctfld_module, only: nOrdInp, nSInit, RadInp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NAt, ICharg, IAt(NAt), ITypRad, NSphInp, m, iPrint
real(kind=wp), intent(in) :: C(3,NAt)
real(kind=wp), intent(inout) :: Alpha
real(kind=wp), intent(out) :: XSph(m), YSph(m), ZSph(m), Rad(m)
integer(kind=iwp), intent(out) :: NOrd(m)
integer(kind=iwp) :: I
real(kind=wp), allocatable :: Chg(:)

! Assign GEPOL sphere positions and radii according to solute atoms nature
select case (ITypRad)
  case (1)
    ! United Atom Topological Model (UATM) radii:
    call mma_allocate(Chg,NAt,label='Chg')
    Chg(:) = Zero
    call UATM(u6,ICharg,NAt,NSinit,m,Rad,Alpha,C,IAt,NOrd,Chg,iPrint)
    call mma_deallocate(Chg)
  case (2)
    ! Pauling radii on each atom:
    do I=1,NAt
      NOrd(I) = I
      Rad(I) = Pauling(IAt(I))
    end do
    Alpha = 1.2_wp
    NSinit = NAt
  case (3)
    ! Sphere radii given in the input
    NOrd(1:NSphInp) = NOrdInp(1:NSphInp)
    Rad(1:NSphInp) = RadInp(1:NSphInp)
    Alpha = 1.2_wp
    NSinit = NSphInp
  case default
    write(u6,'(a)') 'Unrecognized radii type !'
    call Abend()
end select
if (((ITypRad == 2) .or. (ITypRad == 3)) .and. (iPrint > 5)) call PrtCav(u6,ITypRad,NSinit,NOrd,Alpha,Rad)

do I=1,NSinit
  XSph(I) = C(1,NOrd(I))
  YSph(I) = C(2,NOrd(I))
  ZSph(I) = C(3,NOrd(I))
end do
Rad(1:NSinit) = Rad(1:NSinit)*Alpha

return

end subroutine FndSph
