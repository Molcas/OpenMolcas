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

subroutine check_Fthaw(iRC)

use OFembed, only: ThrFThaw
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: iRC
#include "warnings.h"
integer(kind=iwp) :: i, iSeed, iter, iter0, Lu
real(kind=wp) :: DEneA, DEneB, E1, E3, EneA, EneB
logical(kind=iwp) :: ok
real(kind=wp), allocatable :: Ene(:,:)
integer(kind=iwp), external :: IsFreeUnit

if (ThrFThaw <= Zero) return

call f_inquire('AUXRFIL',ok)
if (.not. ok) return

call NameRun('AUXRFIL')
call Get_dScalar('Last energy',EneB)
call NameRun('#Pop')
call Get_dScalar('Last energy',EneA)

iSeed = 7
Lu = IsFreeUnit(iSeed)
call f_inquire('FRETHAW',ok)
if (.not. ok) then
  call molcas_open(Lu,'FRETHAW')
  !open(Lu,file='FRETHAW')
  write(Lu,'(I4,2F18.10)') 1,EneA,EneB
  close(Lu,status='keep')
else
  call molcas_open(Lu,'FRETHAW')
  !open(Lu,file='FRETHAW',status='old')

  read(Lu,'(I4,2F18.10)') iter0,E1,E3
  if (iter0 == 1000) then
    write(u6,*) ' Error! check_Fthaw: maxIter reached! '
    call Abend()
  end if
  call mma_allocate(Ene,iter0,4,label='Ene')
  Ene(1,1) = E1
  Ene(1,3) = E3
  do i=2,iter0
    read(Lu,'(I4,4F18.10)') iter,Ene(i,1),Ene(i,2),Ene(i,3),Ene(i,4)
  end do

  iter = iter0+1
  DEneA = EneA-Ene(iter0,1)
  DEneB = EneB-Ene(iter0,3)

  rewind Lu
  write(Lu,'(I4,2F18.10)') iter,Ene(1,1),Ene(1,3)
  do i=2,iter0
    write(Lu,'(I4,4F18.10)') iter,Ene(i,1),Ene(i,2),Ene(i,3),Ene(i,4)
  end do
  write(Lu,'(I4,4F18.10)') iter,EneA,DEneA,EneB,DEneB

  write(u6,*)
  write(u6,*) '*******************************************************************************'
  write(u6,*) '*************** Energy Statistics for Freeze-n-Thaw ***************************'
  write(u6,*) '*******************************************************************************'
  write(u6,*) '         Energy_A       Delta(Energy_A)      Energy_B       Delta(Energy_B)'
  write(u6,'(I3,1X,F18.10,18X,F18.10)') 1,Ene(1,1),Ene(1,3)
  do i=2,iter0
    write(u6,'(I3,1X,4F18.10)') i,Ene(i,1),Ene(i,2),Ene(i,3),Ene(i,4)
  end do
  write(u6,'(I3,1X,4F18.10)') iter,EneA,DEneA,EneB,DEneB
  write(u6,*) '*******************************************************************************'

  if ((abs(DEneA) < ThrFThaw) .and. (abs(DEneB) < ThrFThaw)) then
    write(u6,'(A,E9.2,A)') ' Convergence reached ! (Thr = ',ThrFThaw,')'
    write(u6,*)
    iRC = _RC_ALL_IS_WELL_
    close(Lu,status='delete')
  else
    write(u6,'(A,E9.2,A)') ' Convergence NOT reached yet ! (Thr = ',ThrFThaw,')'
    write(u6,*)
    close(Lu,status='keep')
  end if
  call mma_deallocate(Ene)
end if

return

end subroutine check_Fthaw
