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

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "rctfld.fh"
parameter(ISAX=1000)
dimension C(3,*), IAt(*)
dimension XSph(*), YSph(*), ZSph(*), Rad(*), NOrd(*)
data IOut/6/,Zero/0.0d0/

! Assign GEPOL sphere positions and radii according to solute atoms nature
if (ITypRad == 1) then
  ! United Atom Topological Model (UATM) radii:
  call GetMem('Chg','Allo','Real',ip_Chg,NAt)
  call dcopy_(NAt,[Zero],0,Work(ip_chg),1)
  call UATM(IOut,ICharg,NAt,NSinit,ToAng,Rad,Alpha,C,IAt,NOrd,Work(ip_Chg),iPrint)
  call GetMem('Chg','Free','Real',ip_Chg,NAt)
elseif (ITypRad == 2) then
  ! Pauling radii on each atom:
  do I=1,NAt
    NOrd(I) = I
    Rad(I) = Pauling(IAt(I))
    Alpha = 1.2
    NSinit = NAt
  end do
elseif (ITypRad == 3) then
  ! Sphere radii given in the input
  do I=1,NSphInp
    NOrd(I) = NOrdInp(I)
    Rad(I) = RadInp(I)
    Alpha = 1.2
    NSinit = NSphInp
  end do
else
  write(6,'(a)') 'Unrecognized radii type !'
  call Abend()
end if
if (((ITypRad == 2) .or. (ITypRad == 3)) .and. (iPrint > 5)) call PrtCav(IOut,ITypRad,NSinit,NOrd,Alpha,Rad)

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

implicit real*8(A-H,O-Z)
dimension NOrd(*), Rad(*)

write(iOut,*)
write(iOut,*)
write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
write(iOut,'(6X,A)') '================================'
if (ITyp == 2) write(iOut,'(6X,A)') 'Pauling radii'
if (ITyp == 3) write(iOut,'(6X,A)') 'Sphere radii from input'
write(iOut,*)
write(IOut,'(6X,A)') ' Nord  Alpha  Radius'
do IS=1,NS
  write(IOut,'(6X,1X,I3,3X,F4.2,3X,F5.3)') NOrd(IS),Alpha,Rad(IS)
end do
write(IOut,'(6X,1X,78("-"))')
write(IOut,*)

return

end subroutine PrtCav
