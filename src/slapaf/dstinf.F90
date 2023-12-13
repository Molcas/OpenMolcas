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

subroutine DstInf(iStop,Just_Frequencies)

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: AtomLbl, Coor, Cx, Dmp_Slapaf, dqInt, Energy, iOptC, iter, lOld_Implicit, Max_Center, MaxItr, MF, mTROld, &
                       Numerical, qInt, RtRnc, SlStop, Weights
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iStop
logical(kind=iwp), intent(in) :: Just_Frequencies
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: i, iDo_dDipM, iIrrep, iOff, iPrint, iRout, isAtom, iTemp, j, jsAtom, LOut, Lu_xyz, N_ZMAT, nCoord, nsAtom_p, &
                     nTemp
real(kind=wp) :: r, r_Iter, x1, x2, xWeight, y1, y2, z1, z2
logical(kind=iwp) :: do_fullprintcoords, do_printcoords, Found
character(len=LenIn), allocatable :: LblTMP(:)
character(len=2), allocatable :: Element(:)
real(kind=wp), allocatable :: CC(:,:), Cx_p(:,:), RV(:,:), xyz(:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 53
iPrint = nPrint(iRout)
do_printcoords = iPrint >= 5
do_fullprintcoords = ((iPrint > 5) .or. (iStop > 1))
LOut = u6
!                                                                      *
!***********************************************************************
!                                                                      *
! Write information of this iteration to the RLXITR file

call Dmp_Slapaf(SlStop,Just_Frequencies,Energy(1),iter,MaxItr,mTROld,lOld_Implicit,size(Coor,2))

if (SuperName /= 'numerical_gradient') then
  call Put_dArray('qInt',qInt,size(qInt))
  call Put_dArray('dqInt',dqInt,size(dqInt))
end if

if (Just_Frequencies) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Geometry information

if (SlStop .or. do_printcoords) then
  write(LOut,*)
  call CollapseOutput(1,'Geometry section')
  write(LOut,*)
  write(LOut,'(80A)') ('*',i=1,80)
  if (SlStop) then
    write(LOut,*) ' Geometrical information of the final structure'
    r_Iter = real(Iter,kind=wp)
    call Add_Info('GEO_ITER',[r_Iter],1,8)
  else if (do_printcoords) then
    write(LOut,*) ' Geometrical information of the new structure'
  end if
  write(LOut,'(80A)') ('*',i=1,80)
  write(LOut,*)
end if

call Get_iScalar('Pseudo atoms',nsAtom_p)
if (nsAtom_p > 0) then
  call mma_allocate(Cx_p,3,nsAtom_p,Label='Cx_p')
  call Get_dArray('Pseudo Coordinates',Cx_p,3*nsAtom_p)
end if

call mma_Allocate(CC,3,nIrrep*(size(Coor,2)+nsAtom_p),Label='CC')
nTemp = 0
call mma_Allocate(LblTMP,nIrrep*(size(Coor,2)+nsAtom_p),Label='LblTMP')
do isAtom=1,size(Coor,2)+nsAtom_p
  if (isAtom <= size(Coor,2)) then
    x1 = Coor(1,isAtom)
    y1 = Coor(2,isAtom)
    z1 = Coor(3,isAtom)
  else
    jsAtom = isAtom-size(Coor,2)
    x1 = Cx_p(1,jsAtom)
    y1 = Cx_p(2,jsAtom)
    z1 = Cx_p(3,jsAtom)
  end if
  doirrep: do iIrrep=0,nIrrep-1
    x2 = x1
    if (btest(iOper(iIrrep),0)) x2 = -x2
    y2 = y1
    if (btest(iOper(iIrrep),1)) y2 = -y2
    z2 = z1
    if (btest(iOper(iIrrep),2)) z2 = -z2

    ! Check if it is already in the list

    do iTemp=1,nTemp
      r = (x2-CC(1,iTemp))**2+(y2-CC(2,iTemp))**2+(z2-CC(3,iTemp))**2
      if (r == Zero) cycle doirrep
    end do
    nTemp = nTemp+1
    if (nTemp > nIrrep*(size(Coor,2)+nsAtom_p)) then
      call WarningMessage(2,'Error in DstInf')
      write(u6,*) 'nTemp > nIrrep*size(Coor,2)'
      call Abend()
    end if
    CC(1,nTemp) = x2
    CC(2,nTemp) = y2
    CC(3,nTemp) = z2
    if (isAtom <= size(Coor,2)) then
      LblTMP(nTemp) = AtomLbl(isAtom)
    else
      LblTMP(nTemp) = 'PC'
    end if
  end do doirrep
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the new cartesian symmetry coordinates.

if (SlStop) then
  write(LOut,*) ' NOTE: on convergence the final predicted structure will be printed here.'
  write(LOut,*) ' This is not identical to the structure printed in the head of the output.'
  call OutCoor('* Nuclear coordinates of the final structure / Bohr     *',AtomLbl,size(Coor,2),Coor,3,size(Coor,2),.false.)
  call OutCoor('* Nuclear coordinates of the final structure / Angstrom *',AtomLbl,size(Coor,2),Coor,3,size(Coor,2),.true.)
else if (Do_PrintCoords) then
  call OutCoor('* Nuclear coordinates for the next iteration / Bohr     *',AtomLbl,size(Coor,2),Coor,3,size(Coor,2),.false.)
  call OutCoor('* Nuclear coordinates for the next iteration / Angstrom *',AtomLbl,size(Coor,2),Coor,3,size(Coor,2),.true.)
end if

if (nsAtom_p > 0) then
  iOff = nTemp-nsAtom_p+1
  call OutCoor('* Pseudo charge coordinates for the next iteration / Bohr     *',LblTMP(iOff),nsAtom_p,Cx_p,3,nsAtom_p,.false.)
  call OutCoor('* Pseudo Charge coordinates for the next iteration / Angstrom *',LblTMP(iOff),nsAtom_p,Cx_p,3,nsAtom_p,.true.)
  call mma_deallocate(Cx_p)
end if

if (do_printcoords) then
  call Get_iScalar('N ZMAT',N_ZMAT)
  if (N_ZMAT > 0) call OutZMAT(size(Coor,2),Coor,N_ZMAT)

  if (do_fullprintcoords) then
    if (nTemp >= 2) call Dstncs(LblTMP,CC,nTemp,Angstrom,Max_Center,5)

    if (nTemp >= 3) call Angles(LblTMP,CC,nTemp,Rtrnc,Max_Center)

    if (nTemp >= 4) call Dihedr(LblTMP,CC,nTemp,Rtrnc,Max_Center)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(LblTMP)
call mma_deallocate(CC)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the new cartesian symmetry coordinates on GEONEW

call Put_Coord_New(Cx(1,1,iter+1),size(Coor,2))
!                                                                      *
!***********************************************************************
!                                                                      *
! If two runfiles are associated with the calculation update both files.

call f_Inquire('RUNFILE2',Found)
if (Found) then
  call NameRun('RUNFILE2')
  call Put_Coord_New(Cx(1,1,iter+1),size(Coor,2))
  call NameRun('#Pop')
end if

! Update the .Opt.xyz file

if (.not. Numerical) then
  call Get_nAtoms_All(nCoord)
  call mma_allocate(xyz,3,nCoord,Label='xyz')
  call mma_allocate(Element,nCoord,Label='Element')
  call Get_Coord_New_All(xyz,nCoord)
  call Get_Name_All(Element)

  Lu_xyz = IsFreeUnit(11)
  call MOLCAS_Open(Lu_xyz,'XYZ')
  write(Lu_xyz,'(I4)') nCoord
  write(Lu_xyz,*) Energy(1+Iter)
  !write(Lu_xyz,'(A)') 'Coordinates generated by Slapaf'
  do i=1,nCoord
    write(Lu_xyz,'(A2,3F15.8)') Element(i),(Angstrom*xyz(j,i),j=1,3)
  end do
  close(Lu_xyz)
  call mma_deallocate(Element)
  call mma_deallocate(xyz)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! If a transition state optimization put the "reaction" vector
! on the RUNFILE(S)

if ((.not. btest(iOptC,7)) .and. SlStop) then

  call mma_allocate(RV,3,size(Coor,2),Label='RV')
  RV(:,:) = MF(:,:)
  do i=1,size(Coor,2)
    xWeight = Weights(i)
    RV(:,i) = RV(:,i)/xWeight
  end do
  call OutCoor('* The Cartesian Reaction vector                         *',AtomLbl,size(Coor,2),RV,3,size(Coor,2),.true.)

  call f_Inquire('RUNREAC',Found)
  if (Found) then
    call NameRun('RUNREAC')
    call Put_dArray('Reaction Vector',RV,3*size(Coor,2))
    call NameRun('#Pop')
  end if
  call f_Inquire('RUNPROD',Found)
  if (Found) then
    call NameRun('RUNPROD')
    call Put_dArray('Reaction Vector',RV,3*size(Coor,2))
    call NameRun('#Pop')
  end if
  call Put_dArray('Reaction Vector',RV,3*size(Coor,2))
  call mma_deallocate(RV)
  iDo_dDipM = 0
  call GF_on_the_fly(iDo_dDipM)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (SlStop .or. do_printcoords) call CollapseOutput(0,'Geometry section')

return

end subroutine DstInf
