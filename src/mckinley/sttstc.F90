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

subroutine Sttstc()

use McKinley_global, only: CPUStat, nFckAcc, nIntegrals, nMOTrans, nOneel, nScreen, nTotal, nTrans, nTwoDens, nTwoel
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iFld
real(kind=wp) :: Diverse, TotCpu
character(len=*), parameter :: NamFld(nTotal+1) = ['1) Calculation of one electron integrals        :', &
                                                   '2) Calculation of two electron integrals        :', &
                                                   '     a) Decontraction of two electron density   :', &
                                                   '     b) Integral evalution 2nd derivatives      :', &
                                                   '     c) Screening                               :', &
                                                   '     d) Transformation of integrals             :', &
                                                   '     e) Direct Fock matrix generation           :', &
                                                   '     f) Direct MO transformation                :', &
                                                   '3)  Control and input                           :', &
                                                   '   T O T A L                                    :']

! Write out timing informations
write(u6,*)
call CollapseOutput(1,'Statistics and timing')
write(u6,'(3X,A)') '---------------------'
write(u6,*)
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,100) '   Part of the program                              CPU    fraction'
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
if (CPUStat(nTotal) > 0.01_wp) then
  TotCpu = CPUStat(nTotal)
else
  TotCpu = 0.01_wp
end if
!TotCpu = max(0.01_wp,CPUStat(nTotal))
CPUStat(nTwoel) = CPUStat(nIntegrals)+CPUStat(nScreen)+CPUStat(nTrans)+CPUStat(nTwoDens)+CPUStat(nFckAcc)+CPUStat(nMOTrans)
Diverse = CPUStat(nTotal)-CPUStat(nTwoEl)-CPUStat(nOneel)

do iFld=1,nTotal-1
  write(u6,'(2x,A45,2f10.2)') NamFld(iFld),CPUStat(iFld),CPUStat(iFld)/TotCpu
end do
write(u6,'(2x,A45,2f10.2)') NamFld(nTotal),diverse,diverse/TotCpu

write(u6,*)
write(u6,'(2x,A45,2F10.2)') NamFld(nTotal+1),TotCpu
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
call CollapseOutput(0,'Statistics and timing')
write(u6,*)

return

100 format(2x,a)

end subroutine Sttstc
