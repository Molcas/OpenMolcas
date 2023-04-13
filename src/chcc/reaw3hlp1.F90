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

subroutine ReaW3hlp1(Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
! this routine does:
! reconstruct  Ww(a",be',b,i)  for aSGrp>beSGrp
! from (a",be'|b,i) records in V3 file LunName

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Ww(1), Wx(1)
integer(kind=iwp) :: dima, dimbe, dimb, no, LunAux
character(len=8) :: LunName
integer(kind=iwp) :: length

!open(unit=LunAux,file=LunName,form='unformatted')
call Molcas_BinaryOpen_Vanilla(LunAux,LunName)

length = dima*dimbe*dimb*no

! read block (a",be'|b,_i)
call rea_chcc(LunAux,length,Ww(1))
!mp call mv0zero(length,length,Ww(1))

close(LunAux)

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Wx)

end subroutine ReaW3hlp1
