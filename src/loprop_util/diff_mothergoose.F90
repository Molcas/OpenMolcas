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

subroutine Diff_MotherGoose(Diffuse,nAt,nB,MP,nij,EC,ANr,Ttot,Ttot_Inv,lMax,TP,dLimmo,Thrs1,Thrs2,nThrs,iPrint,ThrsMul,LuYou)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Diffuse(3)
integer(kind=iwp), intent(in) :: nAt, nB, nij, ANr(nAt), lMax, nThrs, iPrint, LuYou
real(kind=wp), intent(in) :: MP(nij,*), EC(3,nij), Ttot(nB,nB), Ttot_Inv(nB,nB), TP(nAt), dLimmo(2), Thrs1, Thrs2, ThrsMul
integer(kind=iwp) :: lMaxF
real(kind=wp), allocatable :: Pot_Expo(:), Pot_Fac(:), Pot_Point(:)
logical(kind=iwp), allocatable :: Diffed(:)

! Sag hej till publiken.

write(u6,'(A)') '  Enter Slater charge distribution section.'

! Take different route for the different methods for getting
! diffuse distributions.

call mma_allocate(Pot_Expo,nij*2,label='Pot_Expo')
call mma_allocate(Pot_Point,nij,label='Pot_Point')
call mma_allocate(Pot_Fac,nij*4,label='Pot_Fac')
call mma_allocate(Diffed,nij*2,label='Diffed')

if (Diffuse(2)) then
  write(u6,'(A)') '    ---Run a non-linear fit, (Levenberg-Marquart).'
  write(u6,'(A)') '        Thresholds'
  write(u6,991) '           Delta                   :',Thrs1
  write(u6,991) '           Lambda                  :',Thrs2
  write(u6,991) '           Factor                  :',ThrsMul
  write(u6,992) '           Min. decreasing steps   :',nThrs
  write(u6,'(A)') '        Local limit factors'
  write(u6,993) '           Low:',dLimmo(1),'     High:',dLimmo(2)
  call Diff_Numerical(nAt,nB,MP,nij,EC,ANr,Ttot,Ttot_Inv,lMax,TP,dLimmo,Thrs1,Thrs2,nThrs,iPrint,ThrsMul,Pot_Expo,Pot_Point, &
                      Pot_Fac,Diffed)
else if (Diffuse(3)) then
  write(u6,*)
  write(u6,*) 'Not programmed yet, bitte sehr.'
  call Abend()
end if

! Print, analyze uzw, the result of the diffuse stuff.

call WeGotThis(nAt,nB,MP,nij,EC,lMax,iPrint,Pot_Expo,Pot_Point,Pot_Fac,Diffed)

! Generate file with information for other programs.

lMaxF = 1
call YouGetThis(EC,Pot_Expo,Pot_Point,Pot_Fac,Diffed,MP,lMax,lMaxF,nij,LuYou)

call mma_deallocate(Pot_Expo)
call mma_deallocate(Pot_Point)
call mma_deallocate(Pot_Fac)
call mma_deallocate(Diffed)

return

991 format(A,ES12.5)
992 format(A,I2)
993 format(2(A,F10.5))

end subroutine Diff_MotherGoose
