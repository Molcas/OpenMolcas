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

subroutine Diff_MotherGoose(Diffuse,nAt,nB,ipMP,ipC,nij,ip_EC,ip_ANr,ip_Ttot,ip_Ttot_Inv,lMax,iTP,dLimmo,Thrs1,Thrs2,nThrs,iPrint, &
                            ThrsMul,LuYou)

use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp) :: Diffuse(3)
integer(kind=iwp) :: nAt, nB, ipMP, ipC, nij, ip_EC, ip_ANr, ip_Ttot, ip_Ttot_Inv, lMax, iTP, nThrs, iPrint, LuYou
real(kind=wp) :: dLimmo(2), Thrs1, Thrs2, ThrsMul
integer(kind=iwp) :: lMaxF
real(kind=wp) :: Pot_Expo(nij*2), Pot_Fac(nij*4), Pot_Point(nij)
logical(kind=iwp) :: Diffed(nij*2)
#include "WrkSpc.fh"

! Sag hej till publiken.

write(u6,'(A)') '  Enter Slater charge distribution section.'

! Take different route for the different methods for getting
! diffuse distributions.

if (Diffuse(2)) then
  write(u6,'(A)') '    ---Run a non-linear fit, (Levenberg-Marquart).'
  write(u6,'(A)') '        Thresholds'
  write(u6,991) '           Delta                   :',Thrs1
  write(u6,991) '           Lambda                  :',Thrs2
  write(u6,991) '           Factor                  :',ThrsMul
  write(u6,992) '           Min. decreasing steps   :',nThrs
  write(u6,'(A)') '        Local limit factors'
  write(u6,993) '           Low:',dLimmo(1),'     High:',dLimmo(2)
  call Diff_Numerical(nAt,nB,ipMP,ipC,nij,Work(ip_EC),iWork(ip_ANr),ip_Ttot,ip_Ttot_Inv,lMax,iTP,dLimmo,Thrs1,Thrs2,nThrs,iPrint, &
                      ThrsMul,Pot_Expo,Pot_Point,Pot_Fac,Diffed)
elseif (Diffuse(3)) then
  write(u6,*)
  write(u6,*) 'Not programmed yet, bitte sehr.'
  call Abend()
end if
991 format(A,E12.5)
992 format(A,I2)
993 format(2(A,F10.5))

! Print, analyze uzw, the result of the diffuse stuff.

call WeGotThis(nAt,nB,ipMP,ipC,nij,Work(ip_EC),iWork(ip_ANr),ip_Ttot,ip_Ttot_Inv,lMax,iTP,iPrint,Pot_Expo,Pot_Point,Pot_Fac,Diffed)

! Generate file with information for other programs.

lMaxF = 1
call YouGetThis(nAt,Work(ip_EC),Pot_Expo,Pot_Point,Pot_Fac,Diffed,ipMP,lMax,lMaxF,nij,LuYou)

return

end subroutine Diff_MotherGoose
