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

subroutine Print_qEVec2(nH,EVal,EVec)

implicit real*8(a-h,o-z)
real*8 EVec(nH,nH), EVal(nH*(nH+1)/2)
character(len=14) qLbl(nH)
character(len=14) cLbl
character(len=120) Temp

! Skip Primitive Coords

Lu_UDIC = 91
Temp = 'UDIC'
call molcas_open(Lu_UDIC,Temp)
10 read(Lu_UDIC,'(A)') Temp
call UpCase(Temp)
if (Temp(1:4) == 'VARY') Go To 20
goto 10

! Read Internal Coords Labels

20 do iLines=1,nH
40 read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') Go To 40
  cLbl = ' '
  do j=1,14
    if ((Temp(j:j) == ' ') .or. (Temp(j:j) == '=')) goto 30
    cLbl(j:j) = Temp(j:j)
  end do
30 continue
  qLbl(iLines) = cLbl
end do

Lu = 6
IncQQ = 5
do iiQQ=1,nH,IncQQ
  mQQ = min(nH,iiQQ+IncQQ-1)
  write(Lu,*)
  write(Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
  write(Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',(EVal(iQQ*(iQQ+1)/2),iQQ=iiQQ,mQQ)
  write(Lu,*)
  do iq=1,nH
    write(Lu,'(1X,A,5F10.6)') qLbl(iq),(EVec(iq,iQQ),iQQ=iiQQ,mQQ)
  end do
  write(Lu,*)
end do

close(Lu_UDIC)

return

end subroutine Print_qEVec2
