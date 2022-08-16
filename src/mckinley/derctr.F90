!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine DerCtr(mdci,mdcj,mdck,mdcl,ldot,JfGrd,IndGrd,JfHss,IndHss,JfG,mbatch)

use Center_Info
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "disp2.fh"
logical JfHss(4,3,4,3), IfHss(4,3,4,3), JfGrd(3,4), IfGrd(3,4), IfG(4), JfG(4), ldot
integer IndHss(4,3,4,3,0:7), JndHss(4,3,4,3,0:7), IndGrd(3,4,0:7), JndGrd(3,4,0:7)
logical, external :: TF
!define _OLD_CODE_
#ifdef _OLD_CODE_
integer iCo(4), iCom(0:7,0:7), iStabM(0:7), idcrr(0:7)
logical chck
logical, external :: TstFnc
#endif
! Statement function
Ind(i1,i2) = i1*(i1-1)/2+i2

nnIrrep = nIrrep
call lCopy(12,[.false.],0,ifgrd,1)
if (sIrrep) nnIrrep = 1

do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdci,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdci,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,1,iIrrep) = nDisp
      IfGrd(iCar+1,1) = .true.
    else
      IndGrd(iCar+1,1,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdcj,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdcj,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,2,iIrrep) = nDisp
      IfGrd(iCar+1,2) = .true.
    else
      IndGrd(iCar+1,2,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdck,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdck,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,3,iIrrep) = nDisp
      IfGrd(iCar+1,3) = .true.
    else
      IndGrd(iCar+1,3,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdcl,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdcl,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,4,iIrrep) = nDisp
      IfGrd(iCar+1,4) = .true.
    else
      IndGrd(iCar+1,4,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  do iCar=1,3
    do iSh=1,4
      JndGrd(iCar,iSh,iIrrep) = IndGrd(iCar,iSh,iIrrep)
      JfGrd(iCar,iSh) = IfGrd(iCar,iSh)
    end do
  end do
end do
#ifdef _OLD_CODE_
iCo(1) = mdci
iCo(2) = mdcj
iCo(3) = mdck
iCo(4) = mdcl
#endif
call iCopy(144*nirrep,[0],0,IndHss,1)
call iCopy(144*nirrep,[0],0,jndHss,1)
call lCopy(144,[.false.],0,IfHss,1)
call lCopy(144,[.false.],0,JfHss,1)
if (.not. ldot) return

do iAtom=1,4
  do jAtom=1,iAtom

    ! This segment of the code is not really needed.
    ! If turned on it should not do much of a difference.

#   ifdef _OLD_CODE_
    call DCR(LmbdR,dc(iCo(iAtom))%iStab,dc(iCo(iAtom))%nStab,dc(iCo(jAtom))%iStab,dc(iCo(jAtom))%nStab,iDCRR,nDCRR)

    ! Find the stabilizer for A and B

    call Inter(dc(iCo(iAtom))%iStab,dc(iCo(iAtom))%nStab,dc(iCo(jAtom))%iStab,dc(iCo(jAtom))%nStab,iStabM,nStabM)

    ! Generate all possible (left) CoSet
    ! To the stabilizer of A and B

    do iIrrep=0,nIrrep-1
      do jOper=0,nStabM-1
        iCoM(iIrrep,jOper) = ieor(iOper(iIrrep),iStabM(jOper))
      end do
    end do

    ! Order the Coset so we will have the unique ones first
    ! Check uniqueness

    nMax = 1
    do j=1,nIrrep-1
      do i=0,nMax-1
        do ielem=0,nStabM-1
          if (iCoM(i,1) == iCoM(j,ielem)) Go To 435
        end do
      end do

      ! Move unique CoSet

      nMax = nMax+1
      do ielem=0,nStabM-1
        iTmp = iCoM(nMax-1,ielem)
        iCoM(nMax-1,ielem) = iCoM(j,ielem)
        iCoM(j,ielem) = iTmp
      end do
      if (nMax == nIrrep/nStabM) Go To 439
435   continue
    end do
439 continue

    ! Check if the derivative is needed in the present symmetry

    nCoM = nIrrep/nStabM

    do iCar=1,3
      if (iAtom == jAtom) then
        istop = iCar
      else
        iStop = 3
      end if
      do jCar=1,istop
        iComp = ieor(2**(iCar-1),2**(jCar-1))
        Chck = TstFnc(iCoM,0,iComp,nStabM)
        if (Chck) IfHss(iAtom,iCar,jAtom,jCar) = .true.
      end do
    end do
#   endif

   ! Calculate the index for the derivative

    do iIrrep=0,nnIrrep-1
      do iCar=1,3
        if (iAtom == jAtom) then
          istop = iCar
        else
          iStop = 3
        end if
        do jCar=1,istop
          IfHss(iAtom,iCar,jAtom,jCar) = .true.
          if ((jndGrd(iCar,iAtom,iIrrep) > 0) .and. (jndGrd(jCar,jAtom,iIrrep) > 0)) then
            IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = Ind(max(JndGrd(iCar,iAtom,iIrrep),jndGrd(jCar,jAtom,iIrrep)), &
                                                       min(jndGrd(iCar,iAtom,iIrrep),jndGrd(jCar,jAtom,iIrrep)))
          else
            IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0
          end if
        end do
      end do
    end do
  end do
end do

! Scramble the control array for the hessian

call LCopy(4,[.true.],0,Ifg,1)
do iAtom=1,4
  JfG(iAtom) = IfG(iAtom)
  do jAtom=1,iAtom
    do iCar=1,3
      if (iAtom == jAtom) then
        iStop = iCar
      else
        iStop = 3
      end if
      do jCar=1,iStop
        if (iAtom >= jAtom) then
          JfHss(iAtom,iCar,jAtom,jCar) = IfHss(iAtom,iCar,jAtom,jCar)
        else if (iAtom < jAtom) then
          JfHss(iAtom,iCar,jAtom,jCar) = IfHss(jAtom,jCar,iAtom,iCar)
        end if
      end do
    end do
  end do
end do
if (sIrrep) then
  do ii=1,4
    do ic1=1,3
      if (indgrd(ic1,ii,0) == 0) jfgrd(ic1,ii) = .false.
      do ij=1,4
        do ic2=1,3
          if (Indhss(ii,ic1,ij,ic2,0) == 0) JfHss(ii,ic1,ij,ic2) = .false.
        end do
      end do
    end do
  end do
end if
!----------------------------------------------------------------------*
!
!  End Hess
!
!----------------------------------------------------------------------*

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(mbatch)

end subroutine DerCtr
