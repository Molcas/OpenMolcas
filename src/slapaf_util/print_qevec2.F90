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

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(in) :: EVal(nTri_Elem(nH)), EVec(nH,nH)
integer(kind=iwp) :: iiQQ, iLines, IncQQ, iq, iQQ, j, Lu, Lu_UDIC, mQQ
character(len=120) :: Temp
character(len=14) :: cLbl
character(len=14), allocatable :: qLbl(:)
integer(kind=iwp), external :: IsFreeUnit

! Skip Primitive Coords

Lu_UDIC = IsFreeUnit(91)
Temp = 'UDIC'
call molcas_open(Lu_UDIC,Temp)
do
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:4) == 'VARY') exit
end do

! Read Internal Coords Labels

call mma_allocate(qLbl,nH,Label='qLbl')

do iLines=1,nH
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') cycle
  cLbl = ' '
  do j=1,14
    if ((Temp(j:j) == ' ') .or. (Temp(j:j) == '=')) exit
    cLbl(j:j) = Temp(j:j)
  end do
  qLbl(iLines) = cLbl
end do

Lu = u6
IncQQ = 5
do iiQQ=1,nH,IncQQ
  mQQ = min(nH,iiQQ+IncQQ-1)
  write(Lu,*)
  write(Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
  write(Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',(EVal(nTri_Elem(iQQ)),iQQ=iiQQ,mQQ)
  write(Lu,*)
  do iq=1,nH
    write(Lu,'(1X,A,5F10.6)') qLbl(iq),EVec(iq,iiQQ:mQQ)
  end do
  write(Lu,*)
end do

call mma_deallocate(qLbl)

close(Lu_UDIC)

return

end subroutine Print_qEVec2
