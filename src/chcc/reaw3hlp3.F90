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

implicit none
integer dima, dimbe, dimb, no, LunAux
character*8 LunName
real*8 Ww(1:dima,1:dimbe,1:dimb,1:no)
real*8 Wx(*)
! help variables
integer i, a, be, b, beabi, length

! read block (be",a"|b"_i)

length = dima*dimbe*dimb*no
!open(unit=LunAux,file=LunName,form='unformatted')
call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
call rea_chcc(LunAux,length,Wx(1))
!mp call mv0zero(length,length,Wx(1))
close(LunAux)

! Map Ww(a",be",b",i) <- Wx(be",a"|b"_i)

beabi = 0
do i=1,no
  do b=1,dimb
    do a=1,dima
      do be=1,dimbe
        beabi = beabi+1
        Ww(a,be,b,i) = Wx(beabi)
      end do
    end do
  end do
end do

return

end subroutine ReaW3hlp3
