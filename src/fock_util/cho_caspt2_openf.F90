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

subroutine Cho_CASPT2_OpenF(iOpt,iTyp,iSym,nBatch)
! Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
!          delete (iOpt=3) Cholesky vector files for CASPT2 program
!          (full vectors).
!          For iOpt=0, the units are initialized (to -1).
!          iTyp=1: transformed Cholesky vectors (Inactive).
!          iTyp=2: transformed Cholesky vectors (Active).
!          iTyp=3: vectors from (pi|qj) decomposition (Not
!                  implemented yet!)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iOpt, iTyp, iSym, nBatch
#include "WrkSpc.fh"
#include "chocaspt2.fh"
integer(kind=iwp) :: iaddr, iB, LuV, nSym, NCALLS = 0, NUMCHO(8)
character(len=3) :: BaseNm
character(len=7) :: FullNm
character(len=*), parameter :: SecNam = 'Cho_CASPT2_OpenF'

!******************************************************************
if (nBatch > 999) then
  call Cho_x_Quit(SecNam,' nBatch limited to 999 !!!',' ')
end if
call Get_iScalar('nSym',nSym)
call Get_iArray('NumCho',NumCho,nSym)

if (NCALLS == 0) then
  do iB=1,nBatch
    iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
    iWork(iaddr) = -1
  end do
end if

! Initialize units and return for iOpt=0.
! ---------------------------------------

if (iOpt == 0) then
  do iB=1,nBatch
    iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
    iWork(iaddr) = -1
  end do
  return
end if

! Open or close files.
! --------------------
if ((iTyp < 1) .or. (iTyp > 2)) then
  call Cho_x_Quit(SecNam,'iTyp error',' ')
end if

if (iOpt == 1) then
  if (NumCho(iSym) > 0) then
    do iB=1,nBatch
      iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
      if (iWork(iaddr) < 1) then
        call Cho_caspt2_GetBaseNm(BaseNm,iTyp)
        write(FullNm,'(A3,I1,I3)') BaseNm,iSym,iB
        LuV = 7 ! initial guess
        call daName_MF_WA(LuV,FullNm) ! handle inquire/free unit
        iWork(iaddr) = LuV
        write(u6,*) ' Opened file "',FullNm,'" as unit nr LuV=',LuV
        write(u6,*) ' Unit number LuV is stored at address ',iaddr
      end if
    end do
  else
    do iB=1,nBatch
      iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
      iWork(iaddr) = -1
    end do
  end if
else if (iOpt == 2) then
  do iB=1,nBatch
    iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
    if (iWork(iaddr) > 0) then
      write(u6,*) ' Closing lUnit_F=',iWork(iaddr)
      call daClos(iWork(iaddr))
      iWork(iaddr) = -1
    end if
  end do
else if (iOpt == 3) then
  do iB=1,nBatch
    iaddr = ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
    if (iWork(iaddr) > 0) then
      write(u6,*) ' Erasing lUnit_F=',iWork(iaddr)
      call daEras(iWork(iaddr))
      iWork(iaddr) = -1
    end if
  end do
else
  call Cho_x_Quit(SecNam,'iOpt out of bounds',' ')
end if

end subroutine Cho_CASPT2_OpenF
