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

subroutine ReaW3hlp2(Ww,Wx,dima,dimb,no,LunName,LunAux)
! this routine does:
! reconstruct  Ww(a",be',b,i)  for aSGrp=beSGrp
! from (a">=be"|b,i) records in V3 file LunName

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no, LunAux
real(kind=wp), intent(out) :: Ww(dima,dima,dimb,no), Wx(nTri_Elem(dima)*dimb*no)
character(len=8) :: LunName
integer(kind=iwp) :: a, abebi, b, i, length

! read block (a">=be"|b"_i)
length = no*nTri_Elem(dima)*dimb
!open(unit=LunAux,file=LunName,form='unformatted')
call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
call rea_chcc(LunAux,length,Wx)
!mp Wx(1:length) = Zero
close(LunAux)

! Expand and Set Ww(a",be",b",i) <- Wx(a">=be"|b",i)
abebi = 0

do i=1,no
  do b=1,dimb
    do a=1,dima
      Ww(a,1:a-1,b,i) = Wx(abebi+1:abebi+a-1)
      Ww(1:a,a,b,i) = Wx(abebi+1:abebi+a)
      abebi = abebi+a
    end do
  end do
end do

return

end subroutine ReaW3hlp2
