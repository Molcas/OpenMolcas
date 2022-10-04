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

use Symmetry_Info, only: iChTbl, nIrrep, iOper, Symmetry_Info_Get

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
real*8 Vectors(nVectors), Vectors_All(3*nAll_Atoms*nFreq)
integer iGen(3), iCoSet(0:7,0:7), mDisp(0:7), iChCar(3), nDisp(0:7), iStab(0:7)
#ifdef _DEBUGPRINT_
logical Temp
#endif
integer, save :: Active = 0

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
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
write(6,*) 'nGen=',nGen
write(6,*) 'iGen=',(iGen(i),i=1,nGen)
#endif
call ChCar(iChCar,iGen,nGen)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of symmetry adapted displacements

call Get_iScalar('Unique atoms',mUnique_Atoms)
if (mUnique_Atoms /= nUnique_Atoms) then
  write(6,*) 'Get_NMode_All: mUnique_Atoms /= nUnique_Atoms'
  call Abend()
end if
call Allocate_Work(ipCoor,3*mUnique_Atoms)
call Get_dArray('Unique Coordinates',Work(ipCoor),3*mUnique_Atoms)
#ifdef _DEBUGPRINT_
write(6,*) 'nVectors,nAll_Atoms,nFreq=',nVectors,nAll_Atoms,nFreq
write(6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the number of displacements in each irrep

MaxDCR = 0
do iIrrep=0,nIrrep-1
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'iIrrep=',iIrrep
  write(6,*)
# endif
  nDisp(iIrrep) = 0
  ipTmp = ipCoor
  do iUnique_Atom=1,nUnique_Atoms
#   ifdef _DEBUGPRINT_
    write(6,*) 'iUnique_Atom=',iUnique_Atom
#   endif
    iChAtom = iChxyz(Work(ipTmp),iGen,nGen)
    ipTmp = ipTmp+3
    call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
    nCoSet = nIrrep/nStab
#   ifdef _DEBUGPRINT_
    write(6,*) 'nCoSet=',nCoSet
    write(6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
    write(6,*) 'iChAtom=',iChAtom
#   endif
    do iCar=0,2
      iComp = 2**iCar
#     ifdef _DEBUGPRINT_
      write(6,*) 'iComp=',iComp
      Temp = TF(iIrrep,iComp)
      write(6,*) 'TF(iIrrep,iComp)=',Temp
#     endif
      if (TF(iIrrep,iComp)) nDisp(iIrrep) = nDisp(iIrrep)+1
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'Grand Total'
write(6,*) 'nDisp=',(nDisp(i),i=0,nIrrep-1)
write(6,*)
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
  write(6,*) 'iIrrep,nDisp(iIrrep)=',iIrrep,nDisp(iIrrep)
# endif

  do iMode=1,mDisp(iIrrep)
    iFreq = iFreq+1
#   ifdef _DEBUGPRINT_
    write(6,*) 'iMode=',iMode
#   endif

    ! Loop over symmetry unique centers

    ipTmp = ipCoor
    do iUnique_Atom=1,nUnique_Atoms

#     ifdef _DEBUGPRINT_
      write(6,*) 'iUnique_Atom=',iUnique_Atom
#     endif
      ! Get permutational character of the center

#     ifdef _DEBUGPRINT_
      write(6,*) Work(ipTmp),Work(ipTmp+1),Work(ipTmp+2)
#     endif
      iChAtom = iChxyz(Work(ipTmp),iGen,nGen)
      ipTmp = ipTmp+3
      call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
      nCoSet = nIrrep/nStab
#     ifdef _DEBUGPRINT_
      write(6,*) 'nCoSet=',nCoSet
      write(6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
      write(6,*) 'iChAtom=',iChAtom
#     endif

      iVec = 0 ! dummy initialize
      do iCo=0,nCoSet-1
        kOp = iCoSet(iCo,0)
#       ifdef _DEBUGPRINT_
        write(6,*) 'iVector_All=',iVector_All
        write(6,*) 'iVector=',iVector
        write(6,*) 'iCo,kOp=',iCo,kOp
#       endif
        iVec = 0
        do iCar=0,2
          iComp = 2**iCar
          iVector_All = iVector_All+1
#         ifdef _DEBUGPRINT_
          write(6,*) 'iCar=',iCar
          write(6,*) 'iVector_All=',iVector_All
          write(6,*) 'iComp=',iComp
          Temp = TF(iIrrep,iComp)
          !write(6,*) 'TF(iIrrep,iComp)=',Temp
#         endif
          if (TF(iIrrep,iComp)) then
            !write(6,*) 'Belong!'
            iVec = iVec+1

            ! In some cases we only want the normal modes of the first irrep!
            ! In that case branch out.

            if (iVector+iVec > nVectors) exit outer
            Vec = Vectors(iVector+iVec)
            XR = dble(iPrmt(NrOpr(kOp),iComp))
            XY = dble(iChTbl(iIrrep,NrOpr(kOp)))
            Vectors_All(iVector_All) = Vec*XR*XY
          else
            !write(6,*) 'Doesn''t belong!'
            Vectors_All(iVector_All) = Zero
          end if
#         ifdef _DEBUGPRINT_
          write(6,*) 'iVec=',iVec
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
call Free_Work(ipCoor)
#ifdef _DEBUGPRINT_
call RecPrt('Normal mode',' ',Vectors_All,3*nAll_Atoms,nFreq)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

logical function TF(iIrrep,iComp)
  implicit real*8(a-h,o-z)
  logical, external :: TstFnc
  TF = TstFnc(iCoSet,iIrrep,iComp,nIrrep/nCoSet)
end function TF

end subroutine Get_NMode_All
