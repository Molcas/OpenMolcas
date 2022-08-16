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

implicit real*8(a-h,o-z)
#include "cputime.fh"
character*50 NamFld(nTotal+1)
character*60 Fmt
data NamFld /'1) Calculation of one electron integrals        :','2) Calculation of two electron integrals        :', &
             '     a) Decontraction of two electron density   :','     b) Integral evalution  2nd derivatives    :', &
             '     c) Screening                               :','     d) Transfromation of integrals             :', &
             '     e) Direct Fock matrix generation           :','     f) Direct MO transformation                :', &
             '3)  Control and input                           :','   T O T A L                                    :'/

! Write out timing informations
Fmt = '(2x,A)'
write(6,*)
call CollapseOutput(1,'Statistics and timing')
write(6,'(3X,A)') '---------------------'
write(6,*)
write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(6,Fmt) '   Part of the program                              CPU    fraction'
write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
if (CPUStat(nTotal) > 0.01d0) then
  TotCpu = CPUStat(nTotal)
else
  TotCpu = 0.01d0
end if
!TotCpu = max(0.01,CPUStat(nTotal))
CPUStat(nTwoel) = CPUStat(nIntegrals)+CPUStat(nScreen)+CPUStat(nTrans)+CPUStat(nTwoDens)+CPUStat(nFckAck)+CPUStat(nMOTrans)
Diverse = CPUStat(nTotal)-CPUStat(nTwoEl)-CPUStat(nOneel)

do iFld=1,nTotal-1
  write(6,'(2x,A45,2f10.2)') NamFld(iFld),CPUStat(iFld),CPUStat(iFld)/TotCpu
end do
write(6,'(2x,A45,2f10.2)') NamFld(nTotal),diverse,diverse/TotCpu

write(6,*)
write(6,'(2x,A45,2F10.2)') NamFld(nTotal+1),TotCpu
write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
call CollapseOutput(0,'Statistics and timing')
write(6,*)

return

end subroutine Sttstc
