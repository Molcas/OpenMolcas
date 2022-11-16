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

subroutine MkGrid(natom,ipCord,ipGrd,nGrdPt,iRMax,DeltaR,Forces,ipIsMM,iGrdTyp,ipDGrd,nAtQM)

use PCM_arrays, only: Centr, dCntr, DPnt, dRad, dTes, IntSph, NewSph, NVert, PCM_N, PCM_SQ, PCMDM, PCMiSph, PCMSph, PCMTess, SSph, &
                      Vert
use external_centers, only: iXPolType
use stdalloc, only: mma_deallocate
use Constants, only: One, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: natom, ipCord, ipGrd, nGrdPt, iRMax, ipIsMM, iGrdTyp, ipDGrd, nAtQM
real(kind=wp) :: DeltaR
logical(kind=iwp) :: Forces
#include "WrkSpc.fh"
#include "rctfld.fh"
integer(kind=iwp) :: I, iatom, ibla, iCur, ip_LcANr, ip_LcCoor, ipAN, ipBla, ipChrg, ipKeep, iPL, iPnt, iPrint, iPt, ipTmp, &
                     ipTmpDG, ipTmpG, J, jPnt, nDer, New_nGrdPt, nGrdPt_old, nTmp
real(kind=wp) :: R, X, Y, Z
logical(kind=iwp) :: Dirty, Process
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

iPrint = 5
if (iPL >= 3) iPrint = 50
if (iPL >= 4) iPrint = 99
call GetMem('Atomic Numbers','Allo','Inte',ipAN,natom)
call GetMem('Get_Atoms','Allo','Real',ipChrg,natom)
call Get_dArray('Nuclear charge',Work(ipChrg),natom)
do iatom=0,natom-1
  iWork(ipAN+iatom) = int(Work(ipChrg+iatom))
end do
call GetMem('Get_Atoms','Free','Real',ipChrg,natom)
nGrdPt_old = nGrdPt

! PNT grid (it uses Angstroms !!!)

if (abs(iGrdTyp) == 1) then
  DeltaR = DeltaR*Angstrom
  Process = (iGrdTyp == 1)
  call DScal_(3*natom,Angstrom,Work(ipCord),1)
  call PNT(u6,natom,Work(ipCord),iRMax,DeltaR,iWork(ipAN),nGrdPt,Work(ipGrd),iWork(ipIsMM),Process)
  call DScal_(3*natom,One/Angstrom,Work(ipCord),1)
  DeltaR = DeltaR/Angstrom
  if (nGrdPt <= 0) then
    write(u6,'(A)') ' Error in espf/mkgrid: nGrdPt < 0 !!!'
    call Quit_OnUserError()
  end if

  ! Printing the PNT point coordinates

  if (Process .and. (.not. DoDeriv)) then
    if (iPL >= 4) then
      write(u6,'(A,I5,A)') ' PNT grid (in Angstrom) '
      do iPt=1,nGrdPt
        iCur = 3*(iPt-1)
        write(u6,'(A4,3F15.6)') ' X  ',Work(ipGrd+iCur),Work(ipGrd+iCur+1),Work(ipGrd+iCur+2)
      end do
    end if
    call DScal_(3*nGrdPt,One/Angstrom,Work(ipGrd),1)
  end if

else

  ! GEPOL grid (made of iRMax surfaces)

  iXPolType = 0
  PCM = .true.
  DoDeriv = Forces
  nDer = nAtQM*3
  do J=0,iRMax-1
    if (iPL >= 3) write(u6,'(A13,I1)') ' GEPOL shell ',J+1
    call GetMem('LcCoor','Allo','Real',ip_LcCoor,3*natom)
    call GetMem('LcANr','Allo','Inte',ip_LcANr,natom)
    nPCM_info = 0
    call PCM_Cavity(iPrint,0,natom,Work(ipCord),iWork(ipAN),iWork(ipIsMM),Work(ip_LcCoor),iWork(ip_LcANr),J)
    call GetMem('LcANr','Free','Inte',ip_LcANr,natom)
    call GetMem('LcCoor','Free','Real',ip_LcCoor,3*natom)
    if (J == 0) then
      nTmp = 0
      nGrdPt = nTs
      call Allocate_Work(ipTmp,3*nGrdPt)
    else
      nTmp = nGrdPt
      call Allocate_Work(ipBla,3*nTmp)
      call dcopy_(3*nTmp,Work(ipTmp),1,Work(ipBla),1)
      nGrdPt = nTs+nGrdPt
      call Free_Work(ipTmp)
      call Allocate_Work(ipTmp,3*nGrdPt)
      call dcopy_(3*nTmp,Work(ipBla),1,Work(ipTmp),1)
      call Free_Work(ipBla)
    end if
    if (DoDeriv) call GetMem('ESPF_DGrid','Allo','Real',ipDGrd,3*nGrdPt*nDer)
    do I=1,nTs
      call dcopy_(3,PCMTess(1,I),1,Work(ipTmp+3*nTmp+3*(I-1)),1)
    end do
    if (DoDeriv) call dcopy_(3*nGrdPt*nDer,DPnt,1,Work(ipDGrd),1)
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
  call GetMem('ESPF_Grid','Allo','Real',ipGrd,3*nGrdPt)
  call dcopy_(3*nGrdPt,Work(ipTmp),1,Work(ipGrd),1)
  call Free_Work(ipTmp)

  ! Cleaning the GEPOL grid:
  ! all grid points must be distant by 1 bohr at least

  Dirty = .true.
  do while (Dirty)
    Dirty = .false.
    call Allocate_iWork(ipKeep,nGrdPt)
    do iPnt=0,nGrdPt-1
      iWork(ipKeep+iPnt) = 1
    end do
    do iPnt=0,nGrdPt-2
      if (iWork(ipKeep+iPnt) == 0) cycle
      do jPnt=iPnt+1,nGrdPt-1
        X = Work(ipGrd+3*jPnt)-Work(ipGrd+3*iPnt)
        Y = Work(ipGrd+3*jPnt+1)-Work(ipGrd+3*iPnt+1)
        Z = Work(ipGrd+3*jPnt+2)-Work(ipGrd+3*iPnt+2)
        R = sqrt(X*X+Y*Y+Z*Z)
        if (R < One) iWork(ipKeep+jPnt) = 0
      end do
    end do
    New_nGrdPt = 0
    do iPnt=0,nGrdPt-1
      if (iWork(ipKeep+iPnt) == 1) New_nGrdPt = New_nGrdPt+1
    end do
    Dirty = New_nGrdPt < nGrdPt
    if (Dirty) then
      call Allocate_Work(ipTmpG,3*nGrdPt)
      call dcopy_(3*nGrdPt,Work(ipGrd),1,Work(ipTmpG),1)
      call GetMem('ESPF_Grid','Free','Real',ipGrd,3*nGrdPt)
      call GetMem('ESPF_Grid','Allo','Real',ipGrd,3*New_nGrdPt)
      if (DoDeriv) then
        call Allocate_Work(ipTmpDG,3*nGrdPt*NDer)
        call dcopy_(3*nGrdPt*NDer,Work(ipDGrd),1,Work(ipTmpDG),1)
        call GetMem('ESPF_DGrid','Free','Real',ipDGrd,3*nGrdPt*NDer)
        call GetMem('ESPF_DGrid','Allo','Real',ipDGrd,3*New_nGrdPt*NDer)
      end if
      ibla = -1
      do iPnt=0,nGrdPt-1
        if (iWork(ipKeep+iPnt) == 1) then
          ibla = ibla+1
          call dcopy_(3,Work(ipTmpG+3*iPnt),1,Work(ipGrd+3*ibla),1)
          if (DoDeriv) call dcopy_(9*nAtQM,Work(ipTmpDG+iPnt),nGrdPt,Work(ipDGrd+ibla),New_nGrdPt)
        end if
      end do
      if (DoDeriv) call Free_Work(ipTmpDG)
      call Free_Work(ipTmpG)
      nGrdPt = New_nGrdPt
    end if
    call Free_iWork(ipKeep)
  end do

  ! Printing the GEPOL point coordinates

  if ((.not. DoDeriv) .and. (iPL >= 4)) then
    call DScal_(3*nGrdPt,Angstrom,Work(ipGrd),1)
    write(u6,'(A)') 'PCM grid (in Angstroms):'
    do iPnt=0,nGrdPt-1
      iCur = 3*iPnt
      write(u6,'(A4,3F15.6)') ' X  ',Work(ipGrd+iCur),Work(ipGrd+iCur+1),Work(ipGrd+iCur+2)
    end do
    call DScal_(3*nGrdPt,One/Angstrom,Work(ipGrd),1)
  end if
  PCM = .false.
end if
if ((nGrdPt_old /= 0) .and. (nGrdPt /= nGrdPt_old)) then
  write(u6,'(A,2i10)') 'MkGrid: inconsistency in nGrdPt:',nGrdPt_old,nGrdPt
  call Abend()
end if
call GetMem('Atomic Numbers','Free','Inte',ipAN,natom)

return

end subroutine MkGrid
