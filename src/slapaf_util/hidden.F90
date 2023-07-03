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

subroutine Hidden(Coor,AN,nHidden)
! Add to the Grand atom list some hidden atoms, coming e.g.
! from the MM part of a QM/MM system

use Slapaf_Parameters, only: rHidden

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "angstr.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "periodic_table.fh"
real*8, allocatable :: Coor(:,:)
integer, allocatable :: AN(:)
real*8, allocatable :: Coor_h(:,:)
integer, allocatable :: AN_h(:)
logical Do_ESPF, Exist, Exist2
character(len=180) Line
character(len=180), external :: Get_Ln
real*8 XYZ(3)
character(len=2) Symbol
character(len=LenIn), allocatable :: LabMMO(:)
real*8, allocatable :: h_xyz(:,:)
integer, allocatable :: h_AN(:)
!
nHidden = 0
if (rHidden < Two) return
iPL = iPrintLevel(-1)
mTtAtm = size(Coor,2)

!#define _DEBUGPRINT_

#ifdef _DEBUGPRINT_
iPL = 4
#endif
iHidden = 0
Do_ESPF = .false.

! Is there a ESPF/QMMM file ?

call DecideOnESPF(Do_ESPF)
if (Do_ESPF) then

  ! Try MM atoms from Tinker QM/MM interface file

  call F_Inquire('QMMM',Exist)
  if (Exist) then
    ITkQMMM = IsFreeUnit(25)
    call Molcas_Open(ITkQMMM,'QMMM')
    Line = ' '
    do while (index(Line,'TheEnd ') == 0)
      Line = Get_Ln(ITkQMMM)

      ! Read the maximum number of MM atoms + some temporary arrays allocation

      if (index(Line,'NMM') /= 0) then
        call Get_I1(2,nHidden)
        if (iPL > 3) write(6,'(A,I5,A)') 'Found ',nHidden,' hidden atoms.'
        if (nHidden > 0) then
          call mma_allocate(h_xyz,3,nHidden,Label='h_xyz')
          call mma_allocate(h_AN,nHidden,Label='h_AN')
          do iHidden=1,nHidden
            Line = Get_Ln(ITkQMMM)
            if (index(Line,'MMCoord') == 0) then
              write(6,*) 'Error in hidden. Last line does not start with MMCoord:'
              write(6,*) Line
              call Quit_onUserError()
            end if
            call Get_I1(2,iAtNum)
            h_AN(iHidden) = -iAtNum
            call Get_F(3,XYZ,3)
            h_xyz(1,iHidden) = XYZ(1)/Angstr
            h_xyz(2,iHidden) = XYZ(2)/Angstr
            h_xyz(3,iHidden) = XYZ(3)/Angstr
          end do
        end if
      end if
    end do
    close(ITkQMMM)
  end if

  ! Try outer MM atoms stored on runfile

  if (.not. Exist) then
    call Qpg_dArray('MMO Coords',Exist2,nHidden)
  else
    Exist2 = .false.
  end if
  if (Exist2) then
    nHidden = nHidden/3
    call mma_allocate(h_xyz,3,nHidden,Label='h_xyz')
    call mma_allocate(h_AN,nHidden,Label='h_AN')
    call mma_allocate(LabMMO,nHidden,Label='LabMMO')
    call Get_dArray('MMO Coords',h_xyz,nHidden*3)
    call Get_cArray('MMO Labels',LabMMO,LenIn*nHidden)
    do iHidden=1,nHidden
      Symbol(1:1) = LabMMO(iHidden)(1:1)
      Symbol(2:2) = LabMMO(iHidden)(2:2)
      if (Symbol(2:2) == '_') Symbol = ' '//Symbol(1:1)
      do i=0,Num_Elem
        if (Ptab(i) == Symbol) then
          h_AN(iHidden) = -i
          exit
        end if
      end do
    end do
    call mma_deallocate(LabMMO)
  end if
end if
if (iPL > 3) call RecPrt('Hidden coord:',' ',h_xyz,3,nHidden)

! Select the hidden atoms to be kept.

nKept = 0
if (nHidden > 0) call Select_hidden(mTtAtm,nHidden,Coor,h_xyz,h_AN,nKept,iPL)

! Copy all the arrays needed by box and nlm

if (nKept > 0) then
  if (iPL > 3) then
    write(6,'(A8,I5,A)') 'Hidden: ',nKept,' atoms are kept.'
  end if
  mTot = mTtAtm+nKept
  call mma_allocate(Coor_h,3,mTot,Label='Coor_h')
  call mma_allocate(AN_h,mTot,Label='AN_h')

  call dCopy_(3*mTtAtm,Coor,1,Coor_h,1)
  call iCopy(mTtAtm,AN,1,AN_h,1)

  ! Copy the kept hidden atom coordinates, atom numbers and masses

  iKept = 0
  do iHidden=1,nHidden
    if (h_AN(iHidden) > 0) then
      iKept = iKept+1
      call dCopy_(3,h_xyz(:,iHidden),1,Coor_h(:,mTtAtm+iKept),1)
      AN_h(mTtAtm+iKept) = h_AN(iHidden)
    end if
  end do
  if (iKept /= nKept) then
    write(6,'(A)') ' Hidden: wrong number of kept hidden atoms.'
    call Quit_OnUserError()
  end if
  call mma_deallocate(h_AN)
  call mma_deallocate(h_xyz)
  call mma_deallocate(Coor)
  call mma_deallocate(AN)

  ! The end

  call mma_allocate(Coor,3,mTot,Label='Coor')
  Coor(:,:) = Coor_h(:,:)
  call mma_deallocate(Coor_h)
  call mma_allocate(AN,mTot,Label='AN')
  AN(:) = AN_h(:)
  call mma_deallocate(AN_h)

  if (iPL > 3) then
    call RecPrt('Hidden: Coor',' ',Coor,3,mTot)
  end if
end if
nHidden = nKept

return

end subroutine Hidden
