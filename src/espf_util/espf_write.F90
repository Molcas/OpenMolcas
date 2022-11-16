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

subroutine espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,Show_espf,Forces, &
                      DoDirect)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: MltOrd, iRMax, iGrdTyp, nGrdPt, ipMltp, nMult, ipIsMM, natom
logical(kind=iwp) :: DoTinker, DoGromacs, lMorok, Show_espf, Forces, DoDirect
#include "WrkSpc.fh"
integer(kind=iwp) :: iAt, iBlaQ, iMlt, iPL, IPotFl, ITkQMMM, j
real(kind=wp) :: DeltaR, EQMMM
logical(kind=iwp) :: Exists
real(kind=wp), allocatable :: Grad(:,:)
integer(kind=iwp), external :: iPL_espf, IsFreeUnit

! Espf data are saved

iPL = iPL_espf()

! Save data in the ESPF.DATA file

IPotFl = 12
IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.DATA')
write(IPotFl,'(A10,I10)') 'MLTORD    ',MltOrd/4
write(IPotFl,'(A10,I10)') 'IRMAX     ',iRMax
write(IPotFl,'(A10,F12.9)') 'DELTAR    ',DeltaR
write(IPotFl,'(A10,I10)') 'GRIDTYPE  ',iGrdTyp
write(IPotFl,'(A10,I10)') 'GRID      ',nGrdPt
if (DoTinker) write(IPotFl,'(A10)') 'TINKER    '
if (DoGromacs) write(IPotFl,'(A10)') 'GROMACS   '
if (lMorok) write(IPotFl,'(A10)') 'LA_MOROK  '
if (DoDirect) write(IPotFl,'(A10)') 'DIRECT    '
if (ipMltp /= ip_Dummy) then
  write(IPotFl,'(A10,I10)') 'MULTIPOLE ',nMult
  iMlt = 0
  if (MltOrd == 1) then
    do iAt=0,natom-1
      if (iWork(ipIsMM+iAt) == 0) then
        write(IPotFl,'(I6,4F15.8)') iAt+1,Work(ipMltp+iMlt),Zero,Zero,Zero
        iMlt = iMlt+1
      end if
    end do
  else
    do iAt=0,natom-1
      if (iWork(ipIsMM+iAt) == 0) then
        write(IPotFl,'(I6,4F15.8)') iAt+1,(Work(ipMltp+iMlt+j),j=0,3)
        iMlt = iMlt+4
      end if
    end do
  end if
end if
write(IPotFl,'(A10)') 'ENDOFESPF '
close(IPotFl)

if (Show_espf .or. (iPL >= 4)) then
  write(u6,'(/,A,/)') ' Informations found in the ESPF data file:'
  write(u6,'(A10,I10)') ' MLTORD   ',MltOrd/4
  write(u6,'(A10,I10)') ' IRMAX    ',iRMax
  write(u6,'(A10,F12.9)') ' DELTAR   ',DeltaR
  write(u6,'(A10,I10)') ' GRIDTYPE ',iGrdTyp
  write(u6,'(A10,I10)') ' GRID     ',nGrdPt
  if (DoTinker) write(u6,'(A10)') ' TINKER   '
  if (DoGromacs) write(u6,'(A10)') ' GROMACS  '
  if (lMorok) write(u6,'(A10)') ' LA_MOROK '
  if (DoDirect) write(u6,'(A10)') ' DIRECT   '
  if (ipMltp /= ip_Dummy) then
    write(u6,'(A10,I10)') ' MULTIPOLE ',nMult
    iMlt = 0
    if (MltOrd == 1) then
      do iAt=0,natom-1
        if (iWork(ipIsMM+iAt) == 0) then
          write(u6,'(I6,4F15.8)') iAt+1,Work(ipMltp+iMlt),Zero,Zero,Zero
          iMlt = iMlt+1
        end if
      end do
    else
      do iAt=0,natom-1
        if (iWork(ipIsMM+iAt) == 0) then
          write(u6,'(I6,4F15.8)') iAt+1,(Work(ipMltp+iMlt+j),j=0,3)
          iMlt = iMlt+4
        end if
      end do
    end if
  end if
  write(u6,'(A10)') ' ENDOFESPF'
end if

! Special case: Tinker is the driver of the QM/MM calculation.
! QM energy + gradient + ESPF multipoles are stored into the QMMM file

call F_Inquire('QMMM',Exists)
if (Exists .and. Forces .and. (.not. DoTinker)) then
  ITkQMMM = IsFreeUnit(15)
  call Molcas_Open(ITkQMMM,'QMMM')
  call Get_dScalar('Last energy',EQMMM)
  write(ITkQMMM,'(F12.7,I5)') EQMMM,MltOrd/4
  call mma_allocate(Grad,3,nAtom,Label='Grad')
  call Get_dArray_chk('GRAD',Grad,3*nAtom)
  do iAt=1,natom
    iBlaQ = ipMltp+MltOrd*(iAt-1)
    write(ITkQMMM,'(7F12.7)') Grad(1:3,iAt),(Work(iBlaQ+J),J=0,MltOrd-1)
  end do
  close(ITkQMMM)
  call mma_deallocate(Grad)
  close(ITkQMMM)
end if

return

end subroutine espf_write
