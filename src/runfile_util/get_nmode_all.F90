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

subroutine Get_NMode_All(Vectors,nVectors,nFreq,nUnique_Atoms,Vectors_All,nAll_Atoms,mDisp)

use Symmetry_Info, only: iChTbl, iOper, nIrrep, Symmetry_Info_Get
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nVectors, nFreq, nUnique_Atoms, nAll_Atoms, mDisp(0:7)
real(kind=wp), intent(in) :: Vectors(nVectors)
real(kind=wp), intent(_OUT_) :: Vectors_All(3*nAll_Atoms*nFreq)
integer(kind=iwp) :: Active = 0, iC, iCar, iChAtom, iChCar(3), iCo, iComp, iCoSet(0:7,0:7), iFreq, iGen(3), iIrrep, iMode, &
                     iStab(0:7), iUnique_Atom, iVec, iVector, iVector_all, kOp, MaxDCR, mUnique_Atoms, nCoSet, nDisp(0:7), nGen, &
                     nStab
real(kind=wp) :: Vec, XR, XY
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
logical(kind=iwp) :: Temp
#endif
real(kind=wp), allocatable :: Coor(:,:)
integer(kind=iwp), external :: iChxyz, iPrmt, NrOpr

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Vectors',' ',Vectors,1,nVectors)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nGen = 0
if (nIrrep == 2) nGen = 1
if (nIrrep == 4) nGen = 2
if (nIrrep == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen == 3) iGen(3) = iOper(4)
#ifdef _DEBUGPRINT_
write(u6,*) 'nGen=',nGen
write(u6,*) 'iGen=',(iGen(i),i=1,nGen)
#endif
call ChCar(iChCar,iGen,nGen)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of symmetry adapted displacements

call Get_iScalar('Unique atoms',mUnique_Atoms)
if (mUnique_Atoms /= nUnique_Atoms) then
  write(u6,*) 'Get_NMode_All: mUnique_Atoms /= nUnique_Atoms'
  call Abend()
end if
call mma_allocate(Coor,3,mUnique_Atoms,label='Coor')
call Get_dArray('Unique Coordinates',Coor,3*mUnique_Atoms)
#ifdef _DEBUGPRINT_
write(u6,*) 'nVectors,nAll_Atoms,nFreq=',nVectors,nAll_Atoms,nFreq
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the number of displacements in each irrep

MaxDCR = 0
do iIrrep=0,nIrrep-1
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'iIrrep=',iIrrep
  write(u6,*)
# endif
  nDisp(iIrrep) = 0
  iC = 1
  do iUnique_Atom=1,nUnique_Atoms
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iUnique_Atom=',iUnique_Atom
#   endif
    iChAtom = iChxyz(Coor(:,iC),iGen,nGen)
    iC = iC+1
    call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
    nCoSet = nIrrep/nStab
#   ifdef _DEBUGPRINT_
    write(u6,*) 'nCoSet=',nCoSet
    write(u6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
    write(u6,*) 'iChAtom=',iChAtom
#   endif
    do iCar=0,2
      iComp = 2**iCar
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iComp=',iComp
      Temp = TF(iIrrep,iComp)
      write(u6,*) 'TF(iIrrep,iComp)=',Temp
#     endif
      if (TF(iIrrep,iComp)) nDisp(iIrrep) = nDisp(iIrrep)+1
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Grand Total'
write(u6,*) 'nDisp=',(nDisp(i),i=0,nIrrep-1)
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop Irreps

iVector = 0
iVector_all = 0
iFreq = 0
outer: do iIrrep=0,nIrrep-1
# ifdef _DEBUGPRINT_
  write(u6,*) 'iIrrep,nDisp(iIrrep)=',iIrrep,nDisp(iIrrep)
# endif

  do iMode=1,mDisp(iIrrep)
    iFreq = iFreq+1
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iMode=',iMode
#   endif

    ! Loop over symmetry unique centers

    iC = 1
    do iUnique_Atom=1,nUnique_Atoms

#     ifdef _DEBUGPRINT_
      write(u6,*) 'iUnique_Atom=',iUnique_Atom
#     endif
      ! Get permutational character of the center

#     ifdef _DEBUGPRINT_
      write(u6,*) Coor(:,iC)
#     endif
      iChAtom = iChxyz(Coor(:,iC),iGen,nGen)
      iC = iC+1
      call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
      nCoSet = nIrrep/nStab
#     ifdef _DEBUGPRINT_
      write(u6,*) 'nCoSet=',nCoSet
      write(u6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
      write(u6,*) 'iChAtom=',iChAtom
#     endif

      iVec = 0 ! dummy initialize
      do iCo=0,nCoSet-1
        kOp = iCoSet(iCo,0)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'iVector_All=',iVector_All
        write(u6,*) 'iVector=',iVector
        write(u6,*) 'iCo,kOp=',iCo,kOp
#       endif
        iVec = 0
        do iCar=0,2
          iComp = 2**iCar
          iVector_All = iVector_All+1
#         ifdef _DEBUGPRINT_
          write(u6,*) 'iCar=',iCar
          write(u6,*) 'iVector_All=',iVector_All
          write(u6,*) 'iComp=',iComp
          Temp = TF(iIrrep,iComp)
          !write(u6,*) 'TF(iIrrep,iComp)=',Temp
#         endif
          if (TF(iIrrep,iComp)) then
            !write(u6,*) 'Belong!'
            iVec = iVec+1

            ! In some cases we only want the normal modes of the first irrep!
            ! In that case branch out.

            if (iVector+iVec > nVectors) exit outer
            Vec = Vectors(iVector+iVec)
            XR = real(iPrmt(NrOpr(kOp),iComp),kind=wp)
            XY = real(iChTbl(iIrrep,NrOpr(kOp)),kind=wp)
            Vectors_All(iVector_All) = Vec*XR*XY
          else
            !write(u6,*) 'Doesn''t belong!'
            Vectors_All(iVector_All) = Zero
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) 'iVec=',iVec
#         endif
        end do  ! iCar
      end do  ! iCo
      iVector = iVector+iVec

    end do ! iUnique_Atom

#   ifdef _DEBUGPRINT_
    call RecPrt('Normal mode',' ',Vectors_All((iFreq-1)*3*nAll_Atoms+1),3,nAll_Atoms)
#   endif
  end do  ! iMode
end do outer  ! iIrrep
call mma_deallocate(Coor)
#ifdef _DEBUGPRINT_
call RecPrt('Normal mode',' ',Vectors_All,3*nAll_Atoms,nFreq)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

function TF(iIrrep,iComp)

  logical(kind=iwp) :: TF
  integer(kind=iwp), intent(in) :: iIrrep, iComp
  logical(kind=iwp), external :: TstFnc

  TF = TstFnc(iCoSet,iIrrep,iComp,nIrrep/nCoSet)

end function TF

end subroutine Get_NMode_All
