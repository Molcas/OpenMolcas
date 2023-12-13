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
real(kind=wp) :: Coo(3,nCent), CooRef(3,nCent)
integer(kind=iwp) :: i, iCent, iDisk, iLu, ind, indMax, iPart, j, jnd1, jnd2, jP, k, l, nRemoved, nTmp
real(kind=wp) :: dCMx, dCMy, dCMz, dSpread, Esub, Etot, Gamold, GaOld, r, Ract, rMax
logical(kind=iwp) :: Exists, ValidOrNot
character(len=200) :: Head
character(len=6) :: FilSlut, Filstart
real(kind=wp), allocatable :: C(:,:), C2(:,:), Coord(:,:)
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
call mma_allocate(C,nPart*nCent,3,label='Coordinates')
call dDaFile(iLu,2,C,3*nPart*nCent,iDisk)
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
      r = C((j-1)*nCent+1,1)**2+C((j-1)*nCent+1,2)**2+C((j-1)*nCent+1,3)**2
      if (r > rMax) then
        rMax = r
        indMax = j
      end if
    end do
    do j=indMax,nPart-i
      k = (j-1)*nCent
      C(k+1:k+nCent,:) = C(k+nCent+1:k+2*nCent,:)
    end do
  end do

  ! Print the new set of coordinates to a startfile.

  nTmp = nPart-nDel
  iLu = 74
  write(FilSlut,'(A5,i1.1)') 'STFIL',NrStartu
  call DaName(iLu,FilSlut)
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nTmp,Gamold,Gaold,Esub)
  iTcSim(1) = iDisk
  do l=1,3
    call dDaFile(iLu,1,C(:,l),(nPart-nDel)*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nTmp,Gamold,Gaold,Esub)
  call DaClos(iLu)

  ! If user wants, print print print.

  if (iPrint >= 10) then
    do i=1,(nPart-nDel)*nCent
      Coord(:,i) = C(i,:)
    end do
    write(Head,*) 'Final coordinates'
    call Cooout(Head,Coord,nPart-nDel,nCent)
  end if
end if

! If adding solvent molecules.

if (DelOrAdd(2)) then
  do i=1,nPart*nCent
    Coord(:,i) = C(i,:)
  end do
  if (nAdd /= 0) then
    ! Just an ugly trick for using nypart. It requires that the first
    ! slot contains the solvent coordinates, so we, temporarily, put
    ! them there.
    Coord(:,1:nCent) = Cordst(:,1:nCent)
    ! Introduce the new particles. nPart is redefined.
    call NyPart(nAdd,nPart,Coord,rStart,nCent,iSeed)
    ! The ugly trick is reversed, and the first slot is retained.
    do j=1,3
      Coord(j,1:nCent) = C(1:nCent,j)
    end do
  end if

  ! Then dump new coordinates on designated startfile.

  call mma_allocate(C2,nPart*nCent,3,label='NewCoo')
  do i=1,3
    C2(:,i) = Coord(i,:)
  end do
  iLu = 74
  write(FilSlut,'(A5,i1.1)') 'STFIL',NrStartu
  call DaName(iLu,FilSlut)
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,rStart,nPart,Gamold,Gaold,Esub)
  iTcSim(1) = iDisk
  do l=1,3
    call dDaFile(iLu,1,C2(:,l),nPart*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,rStart,nPart,Gamold,Gaold,Esub)
  call DaClos(iLu)
  if (iPrint >= 10) then
    write(Head,*) 'Final coordinates'
    call Cooout(Head,Coord,nPart,nCent)
  end if
  call mma_deallocate(C2)
end if

! If requested, substitute all particles that are not of valid water
! geometry for, you guessed it, valid water molecules.

if (DelOrAdd(3)) then
  nRemoved = 0
  do iPart=1,nPart
    ind = nCent*(iPart-1)
    do iCent=1,nCent
      Coo(:,iCent) = C(ind+iCent,:)
      CooRef(:,iCent) = Cordst(:,iCent)
    end do
    call IsItValid(Coo,CooRef,ValidOrNot)
    if (.not. ValidOrNot) then
      dCMx = Zero
      dCMy = Zero
      dCMz = Zero
      do iCent=1,nCent
        dCMx = dCMx+C(ind+iCent,1)
        dCMy = dCMy+C(ind+iCent,2)
        dCMz = dCMz+C(ind+iCent,3)
      end do
      dCMx = dCMx/real(nCent,kind=wp)
      dCMy = dCMy/real(nCent,kind=wp)
      dCMz = dCMz/real(nCent,kind=wp)
      ! Check if the points are spread out, otherwise just delete.
      dSpread = Zero
      do iCent=1,nCent
        dSpread = dSpread+(C(ind+iCent,1)-dCMx)**2+(C(ind+iCent,2)-dCMy)**2+(C(ind+iCent,3)-dCMz)**2
      end do
      dSpread = dSpread/real(nCent,kind=wp)
      if (dSpread < ThrdSpread) then
        nRemoved = nRemoved+1
        do jP=iPart,nPart-1
          jnd1 = (jP-1)*nCent
          jnd2 = (jP)*nCent
          do iCent=1,nCent
            C(jnd1+iCent,:) = C(jnd2+iCent,:)
          end do
        end do
      else
        do iCent=1,nCent
          C(ind+iCent,1) = dCMx+CooRef(1,iCent)
          C(ind+iCent,2) = dCMy+CooRef(2,iCent)
          C(ind+iCent,3) = dCMz+CooRef(3,iCent)
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
    call dDaFile(iLu,1,C(:,l),nPart*nCent,iDisk)
    iTcSim(1+l) = iDisk
  end do
  iDisk = 0
  call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,Gaold,Esub)
  call DaClos(iLu)

  ! If user wants, print print print.

  if (iPrint >= 10) then
    do i=1,3
      Coord(i,:) = C(:,i)
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
    CooRef(:,1:nCent) = Cordst(:,1:nCent)
    call MoldenDump(C,CooRef,nPart,nCent)
  end if
end if

! Deallocate.

call mma_deallocate(C)

! This routine ends now!

return

end subroutine EditStart
