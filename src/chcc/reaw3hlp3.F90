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

subroutine ReaW3hlp3(Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
! this routine does:
! reconstruct  Ww(a",be',b,i)  for aSGrp<beSGrp
! from (be",a"|b,i) records in V3 file LunName

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, dimb, no, LunAux
real(kind=wp), intent(out) :: Ww(dima,dimbe,dimb,no), Wx(dima*dimbe*dimb*no)
character(len=8) :: LunName
integer(kind=iwp) :: a, b, beabi, i, length

! read block (be",a"|b"_i)

length = dima*dimbe*dimb*no
!open(unit=LunAux,file=LunName,form='unformatted')
call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
call rea_chcc(LunAux,length,Wx)
!mp Wx(1:length) = Zero
close(LunAux)

! Map Ww(a",be",b",i) <- Wx(be",a"|b"_i)

beabi = 0
do i=1,no
  do b=1,dimb
    do a=1,dima
      Ww(a,:,b,i) = Wx(beabi+1:beabi+dimbe)
      beabi = beabi+dimbe
    end do
  end do
end do

return

end subroutine ReaW3hlp3
