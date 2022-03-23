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

subroutine EditStart()

use qmstat_global, only: cDumpForm, Cordst, DelOrAdd, iPrint, iSeed, iTcSim, nAdd, nCent, nDel, nPart, NrStarti, NrStartu, rStart
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "WrkSpc.fh"
integer(kind=iwp) :: iC(3), iC2(3)
real(kind=wp) :: Coo(3,nCent), CooRef(3,nCent)
integer(kind=iwp) :: i, iCent, iDisk, iLu, ind, indMax, iPart, j, jnd1, jnd2, jP, k, l, ll, nRemoved
real(kind=wp) :: dCMx, dCMy, dCMz, dSpread, Esub, Etot, Gamold, GaOld, r, Ract, rMax
logical(kind=iwp) :: Exists, ValidOrNot
character(len=200) :: Head
character(len=6) :: FilSlut, Filstart
real(kind=wp), allocatable :: Coord(:,:)
real(kind=wp), parameter :: ThrdSpread = One
#include "warnings.h"

! Inquire if file exists, and if so open it.

write(FilStart,'(A5,i1.1)') 'STFIL',NrStarti
call f_Inquire(FilStart,Exists)
if (.not. Exists) then
  write(u6,*)
  write(u6,*) 'The input startfile given in the EDITstartfile section was not found.'
  call Quit(_RC_IO_ERROR_READ_)
end if
iLu = 73
call DaName(iLu,Filstart)

! Read header and coordinates.

iDisk = 0
call WrRdSim(iLu,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,GaOld,Esub)
iDisk = iTcSim(1)
do l=1,3
  call GetMem('Coordinates','Allo','Real',iC(l),nPart*nCent)
  call dDaFile(iLu,2,Work(iC(l)),nPart*nCent,iDisk)
end do
call DaClos(iLu)

! Now take different paths depending of what user have requested in the input.

call mma_allocate(Coord,3,nPart*nCent,label='Coord')

! If deleting solvent molecules.

if (DelOrAdd(1)) then

  ! Find the solvent molecules farthest away from origo and delete them.

  do i=1,nDel
    rMax = Zero
    indMax = 0
    do j=1,nPart
      r = Zero
      do k=1,3
        r = r+Work(iC(k)+(j-1)*nCent)**2
      end do
      if (r > rMax) then
        rMax = r
        indMax = j
      end if
    end do
    do j=indMax,nPart-i
      do l=1,3
        do ll=1,nCent
          Work(iC(l)+(j-1)*nCent+ll-1) = Work(iC(l)+j*nCent+ll-1)
        end do
      end do
    end do
  end do

  ! Print the new set of coordinates to a startfile.

  iLu = 74
  write(FilSlut,'(A5,i1.1)') 'STFIL',NrStartu
  call DaName(iLu,FilSlut)
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart-nDel,Gamold,Gaold,Esub)
  iTcSim(1) = iDisk
  do l=1,3
    call dDaFile(iLu,1,Work(iC(l)),(nPart-nDel)*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart-nDel,Gamold,Gaold,Esub)
  call DaClos(iLu)

  ! If user wants, print print print.

  if (iPrint >= 10) then
    do i=1,(nPart-nDel)*nCent
      Coord(1,i) = Work(iC(1)+i-1)
      Coord(2,i) = Work(iC(2)+i-1)
      Coord(3,i) = Work(iC(3)+i-1)
    end do
    write(Head,*) 'Final coordinates'
    call Cooout(Head,Coord,nPart-nDel,nCent)
  end if
end if

! If adding solvent molecules.

if (DelOrAdd(2)) then
  do i=1,nPart*nCent
    Coord(1,i) = Work(iC(1)+i-1)
    Coord(2,i) = Work(iC(2)+i-1)
    Coord(3,i) = Work(iC(3)+i-1)
  end do
  if (nAdd /= 0) then
    ! Just an ugly trick for using nypart. It requires that the first
    ! slot contains the solvent coordinates, so we, temporarily, put
    ! them there.
    do i=1,nCent
      do j=1,3
        Coord(j,i) = Cordst(j,i)
      end do
    end do
    ! Introduce the new particles. nPart is redefined.
    call NyPart(nAdd,nPart,Coord,rStart,nCent,iSeed)
    ! The ugly trick is reversed, and the first slot is retained.
    do i=1,nCent
      do j=1,3
        Coord(j,i) = Work(iC(j)+i-1)
      end do
    end do
  end if

  ! Then dump new coordinates on designated startfile.

  do k=1,3
    call GetMem('NewCoo','Allo','Real',iC2(k),nPart*nCent)
  end do
  do i=1,nPart*nCent
    Work(iC2(1)+i-1) = Coord(1,i)
    Work(iC2(2)+i-1) = Coord(2,i)
    Work(iC2(3)+i-1) = Coord(3,i)
  end do
  iLu = 74
  write(FilSlut,'(A5,i1.1)') 'STFIL',NrStartu
  call DaName(iLu,FilSlut)
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,rStart,nPart,Gamold,Gaold,Esub)
  iTcSim(1) = iDisk
  do l=1,3
    call dDaFile(iLu,1,Work(iC2(l)),nPart*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,rStart,nPart,Gamold,Gaold,Esub)
  call DaClos(iLu)
  if (iPrint >= 10) then
    write(Head,*) 'Final coordinates'
    call Cooout(Head,Coord,nPart,nCent)
  end if
  do k=1,3
    call GetMem('NewCoo','Free','Real',iC2(k),nPart*nCent)
  end do
end if

! If requested, substitute all particles that are not of valid water
! geometry for, you guessed it, valid water molecules.

if (DelOrAdd(3)) then
  nRemoved = 0
  do iPart=1,nPart
    ind = nCent*(iPart-1)
    do iCent=1,nCent
      Coo(1,iCent) = Work(iC(1)+ind+iCent-1)
      Coo(2,iCent) = Work(iC(2)+ind+iCent-1)
      Coo(3,iCent) = Work(iC(3)+ind+iCent-1)
      CooRef(1,iCent) = Cordst(1,iCent)
      CooRef(2,iCent) = Cordst(2,iCent)
      CooRef(3,iCent) = Cordst(3,iCent)
    end do
    call IsItValid(Coo,CooRef,ValidOrNot)
    if (.not. ValidOrNot) then
      dCMx = Zero
      dCMy = Zero
      dCMz = Zero
      do iCent=1,nCent
        dCMx = dCMx+Work(iC(1)+ind+iCent-1)
        dCMy = dCMy+Work(iC(2)+ind+iCent-1)
        dCMz = dCMz+Work(iC(3)+ind+iCent-1)
      end do
      dCMx = dCMx/real(nCent,kind=wp)
      dCMy = dCMy/real(nCent,kind=wp)
      dCMz = dCMz/real(nCent,kind=wp)
      ! Check if the points are spread out, otherwise just delete.
      dSpread = Zero
      do iCent=1,nCent
        dSpread = dSpread+(Work(iC(1)+ind+iCent-1)-dCMx)**2
        dSpread = dSpread+(Work(iC(2)+ind+iCent-1)-dCMy)**2
        dSpread = dSpread+(Work(iC(3)+ind+iCent-1)-dCMz)**2
      end do
      dSpread = dSpread/real(nCent,kind=wp)
      if (dSpread < ThrdSpread) then
        nRemoved = nRemoved+1
        do jP=iPart,nPart-1
          jnd1 = (jP-1)*nCent
          jnd2 = (jP)*nCent
          do iCent=1,nCent
            Work(iC(1)+jnd1+iCent-1) = Work(iC(1)+jnd2+iCent-1)
            Work(iC(2)+jnd1+iCent-1) = Work(iC(2)+jnd2+iCent-1)
            Work(iC(3)+jnd1+iCent-1) = Work(iC(3)+jnd2+iCent-1)
          end do
        end do
      else
        do iCent=1,nCent
          Work(iC(1)+ind+iCent-1) = dCMx+CooRef(1,iCent)
          Work(iC(2)+ind+iCent-1) = dCMy+CooRef(2,iCent)
          Work(iC(3)+ind+iCent-1) = dCMz+CooRef(3,iCent)
        end do
      end if
    end if
  end do
  nPart = nPart-nRemoved

  ! Print the new set of coordinates to a startfile.

  iLu = 74
  write(FilSlut,'(A5,i1.1)') 'STFIL',NrStartu
  call DaName(iLu,FilSlut)
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,Gaold,Esub)
  iTcSim(1) = iDisk
  do l=1,3
    call dDaFile(iLu,1,Work(iC(l)),nPart*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,Gaold,Esub)
  call DaClos(iLu)

  ! If user wants, print print print.

  if (iPrint >= 10) then
    do i=1,nPart*nCent
      Coord(1,i) = Work(iC(1)+i-1)
      Coord(2,i) = Work(iC(2)+i-1)
      Coord(3,i) = Work(iC(3)+i-1)
    end do
    write(Head,*) 'Final coordinates'
    call Cooout(Head,Coord,nPart,nCent)
  end if
end if

call mma_deallocate(Coord)

! If the user want to, print the coordinates in some format suitable
! for graphical representation.

if (DelOrAdd(4)) then
  if (cDumpForm(1:4) == 'MOLD') then
    do iCent=1,nCent
      CooRef(1,iCent) = Cordst(1,iCent)
      CooRef(2,iCent) = Cordst(2,iCent)
      CooRef(3,iCent) = Cordst(3,iCent)
    end do
    call MoldenDump(iC,CooRef,nPart,nCent)
  end if
end if

! Deallocate.

do l=1,3
  call GetMem('Coordinates','Free','Real',iC(l),nPart*nCent)
end do

! This routine ends now!

return

end subroutine EditStart
