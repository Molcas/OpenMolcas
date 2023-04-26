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

subroutine MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,Forces,IsMM,iGrdTyp,DGrid,nAtQM)

use PCM_arrays, only: Centr, dCntr, DPnt, dRad, dTes, IntSph, NewSph, NVert, PCM_N, PCM_SQ, PCMDM, PCMiSph, PCMSph, PCMTess, SSph, &
                      Vert
use external_centers, only: iXPolType
use Data_Structures, only: Alloc2DArray_Type, Alloc4DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom, iRMax, IsMM(natom), iGrdTyp, nAtQM
real(kind=wp), intent(in) :: Cord(3,natom), DeltaR
type(Alloc2DArray_Type), intent(out) :: Grid
integer(kind=iwp), intent(inout) :: nGrdPt
logical(kind=iwp), intent(in) :: Forces
type(Alloc4DArray_Type), intent(out) :: DGrid
#include "rctfld.fh"
integer(kind=iwp) :: ibla, iPL, iPnt, iPrint, iPt, J, jPnt, New_nGrdPt, nGrdPt_old, nTmp
real(kind=wp) :: Dum(1), R, X, Y, Z
logical(kind=iwp) :: Dirty, Process
integer(kind=iwp), allocatable :: AN(:), Keep(:), LcANr(:)
real(kind=wp), allocatable :: Bla(:,:), Chrg(:), LcCoor(:,:), Tmp(:,:), TmpDG(:,:,:,:), TmpG(:,:)
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

iPrint = 5
if (iPL >= 3) iPrint = 50
if (iPL >= 4) iPrint = 99
call mma_allocate(AN,natom,label='Atomic Numbers')
call mma_allocate(Chrg,natom,label='Get_Atoms')
call Get_dArray('Nuclear charge',Chrg,natom)
AN(:) = int(Chrg)
call mma_deallocate(Chrg)
nGrdPt_old = nGrdPt

! PNT grid

if (abs(iGrdTyp) == 1) then
  Process = (iGrdTyp == 1)
  if (Process) then
    call mma_allocate(Grid%A,3,nGrdPt,label='ESPF_Grid')
    call PNT(u6,natom,Cord,iRMax,DeltaR,AN,nGrdPt,Grid%A,IsMM,Process)
  else
    call PNT(u6,natom,Cord,iRMax,DeltaR,AN,nGrdPt,Dum,IsMM,Process)
  end if
  if (nGrdPt <= 0) then
    write(u6,'(A)') ' Error in espf/mkgrid: nGrdPt < 0 !!!'
    call Quit_OnUserError()
  end if

  ! Printing the PNT point coordinates

  if (Process .and. (.not. DoDeriv)) then
    if (iPL >= 4) then
      write(u6,'(A,I5,A)') ' PNT grid (in Angstrom) '
      do iPt=1,nGrdPt
        write(u6,'(A4,3F15.6)') ' X  ',Grid%A(:,iPt)*Angstrom
      end do
    end if
  end if

else

  ! GEPOL grid (made of iRMax surfaces)

  iXPolType = 0
  PCM = .true.
  DoDeriv = Forces
  nGrdPt = 0
  call mma_allocate(Tmp,3,nGrdPt,label='Tmp')
  do J=0,iRMax-1
    if (iPL >= 3) write(u6,'(A13,I1)') ' GEPOL shell ',J+1
    call mma_allocate(LcCoor,3,natom,label='LcCoor')
    call mma_allocate(LcANr,natom,label='LcANr')
    nPCM_info = 0
    call PCM_Cavity(iPrint,0,natom,Cord,AN,IsMM,LcCoor,LcANr,J)
    call mma_deallocate(LcCoor)
    call mma_deallocate(LcANr)
    nTmp = nGrdPt
    nGrdPt = nGrdPt+nTs
    call mma_allocate(Bla,3,nGrdPt,label='Bla')
    Bla(:,1:nTmp) = Tmp
    call move_alloc(Bla,Tmp)
    if (DoDeriv) call mma_allocate(DGrid%A,nGrdPt,nAtQM,3,3,label='ESPF_DGrid')
    Tmp(:,nTmp+1:nTmp+nTs-1) = PCMTess(1:3,:)
    if (DoDeriv) DGrid%A(:,:,:,:) = DPnt
    call mma_deallocate(NewSph)
    call mma_deallocate(IntSph)
    call mma_deallocate(NVert)
    call mma_deallocate(PCMiSph)
    call mma_deallocate(PCM_N)
    call mma_deallocate(PCMDM)
    call mma_deallocate(SSph)
    call mma_deallocate(Centr)
    call mma_deallocate(Vert)
    call mma_deallocate(PCMTess)
    call mma_deallocate(PCMSph)
    if (DoDeriv) then
      call mma_deallocate(dTes)
      call mma_deallocate(dPnt)
      call mma_deallocate(dRad)
      call mma_deallocate(dCntr)
      call mma_deallocate(PCM_SQ)
    end if
  end do
  call mma_allocate(Grid%A,3,nGrdPt,label='ESPF_Grid')
  Grid%A(:,:) = Tmp
  call mma_deallocate(Tmp)

  ! Cleaning the GEPOL grid:
  ! all grid points must be distant by 1 bohr at least

  Dirty = .true.
  do while (Dirty)
    Dirty = .false.
    call mma_allocate(Keep,nGrdPt,label='Keep')
    Keep(:) = 1
    do iPnt=1,nGrdPt-1
      if (Keep(iPnt) == 0) cycle
      do jPnt=iPnt+1,nGrdPt
        X = Grid%A(1,jPnt)-Grid%A(1,iPnt)
        Y = Grid%A(2,jPnt)-Grid%A(2,iPnt)
        Z = Grid%A(3,jPnt)-Grid%A(3,iPnt)
        R = sqrt(X*X+Y*Y+Z*Z)
        if (R < One) Keep(jPnt) = 0
      end do
    end do
    New_nGrdPt = 0
    do iPnt=1,nGrdPt
      if (Keep(iPnt) == 1) New_nGrdPt = New_nGrdPt+1
    end do
    Dirty = New_nGrdPt < nGrdPt
    if (Dirty) then
      call mma_allocate(TmpG,3,New_nGrdPt,label='TmpG')
      if (DoDeriv) call mma_allocate(TmpDG,New_nGrdPt,nAtQM,3,3,label='TmpDG')
      ibla = -1
      do iPnt=1,nGrdPt
        if (Keep(iPnt) == 1) then
          ibla = ibla+1
          Tmp(:,ibla) = Grid%A(:,iPnt)
          if (DoDeriv) TmpDG(ibla,:,:,:) = DGrid%A(iPnt,:,:,:)
        end if
      end do
      call mma_deallocate(Grid%A)
      call move_alloc(TmpG,Grid%A)
      if (DoDeriv) then
        call mma_deallocate(DGrid%A)
        call move_alloc(TmpDG,DGrid%A)
      end if
      nGrdPt = New_nGrdPt
    end if
    call mma_deallocate(Keep)
  end do

  ! Printing the GEPOL point coordinates

  if ((.not. DoDeriv) .and. (iPL >= 4)) then
    write(u6,'(A)') 'PCM grid (in Angstroms):'
    do iPnt=0,nGrdPt-1
      write(u6,'(A4,3F15.6)') ' X  ',Grid%A(:,iPnt)*Angstrom
    end do
  end if
  PCM = .false.
end if
if ((nGrdPt_old /= 0) .and. (nGrdPt /= nGrdPt_old)) then
  write(u6,'(A,2i10)') 'MkGrid: inconsistency in nGrdPt:',nGrdPt_old,nGrdPt
  call Abend()
end if
call mma_deallocate(AN)

return

end subroutine MkGrid
