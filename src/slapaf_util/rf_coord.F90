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

!#define _DEBUGPRINT_
subroutine RF_Coord(nq,nsAtom,iIter,nIter,Cx,Process,Valu,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,Proc_dB,mB_Tot,mdB_Tot,BM,dBM,iBM, &
                    idBM,nB_Tot,ndB_Tot,nqB)

use Symmetry_Info, only: iOper, nIrrep, VarR, VarT
use Slapaf_Info, only: dMass, iCoSet, nStab
use ddvdt, only: Rot_Const, Trans_Const
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, iRef, LuIC, nB_Tot, ndB_Tot
integer(kind=iwp), intent(inout) :: nq, Indq(3,nB), mB_Tot, mdB_Tot, iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter)
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Valu(nB,nIter), fconst(nB), rMult(nB), BM(nB_Tot), dBM(ndB_Tot)
character(len=14), intent(inout) :: qLbl(nB)
#include "print.fh"
integer(kind=iwp) :: i, iAtom, iCent, iDeg, iPrint, iRout, iSym, iTest, ixyz, jxyz, kxyz, mB, nCent, nMass, nOrder, nqRF
real(kind=wp) :: COM_xyz, Deg, RotAng, RotMat(3,3), RotVec(3), TMass, Trans(3), Val
logical(kind=iwp) :: Invariant, PSPrint
character(len=14) :: Label
integer(kind=iwp), allocatable :: iDCR(:), Ind(:)
real(kind=wp), allocatable :: currXYZ(:,:), d2RV(:,:,:), dRVdxyz(:,:,:), Grad(:,:), Hess(:,:), Ref123(:,:), xMass(:)
character(len=*), parameter :: TR_type(6) = ['Tx ','Ty ','Tz ','Ryz','Rzx','Rxy']

iRout = 151
iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
iPrint = 99
#endif

if ((.not. VarR) .and. (.not. VarT)) return
!                                                                      *
!***********************************************************************
!                                                                      *
nqRF = 0
PSPrint = .false.
if (iPrint >= 99) PSPrint = .true.
if (PSPrint) write(u6,*) ' Enter RF_Coords.'

! Find nCent and allocate

nCent = 0
do iAtom=1,nsAtom
  nCent = nCent+nIrrep/nStab(iAtom)
end do
mB = nCent*3
call mma_allocate(currXYZ,3,nCent,label='currXYZ')
call mma_allocate(Ref123,3,nCent,label='Ref123')
call mma_allocate(Grad,3,nCent,label='Grad')
call mma_allocate(dRVdxyz,3,3,nCent,label='dRVdxyz')
call mma_allocate(xMass,nCent,label='xMass')
call mma_allocate(Ind,nCent,label='Ind')
call mma_allocate(iDCR,nCent,label='iDCR')
call mma_allocate(Hess,mB,mB,label='Hess')

! Find index of RF center (origin), etc

iCent = 0
do iAtom=1,nsAtom

  do i=0,nIrrep/nStab(iAtom)-1
    iCent = iCent+1
    call OA(iCoSet(i,iAtom),Cx(:,iAtom,iIter),CurrXYZ(:,iCent))
    call OA(iCoSet(i,iAtom),Cx(:,iAtom,iRef),Ref123(:,iCent))

    Ind(iCent) = iAtom
    iDCR(iCent) = iCoSet(i,iAtom)
  end do
end do

!Fact = One
!if (.not. VarR) Fact = 2.0e-2_wp

!write(u6,*) 'nCent=',nCent
!write(u6,*) (Ind(iCent),iCent=1,nCent)

TMass = Zero
do iCent=1,nCent
  iAtom = Ind(iCent)
  xMass(iCent) = dMass(iAtom)
  TMass = TMass+dMass(iAtom)
end do
! Loop over cartesian components

do ixyz=1,3

  Invariant = .false.
  iTest = 2**(ixyz-1)
  do iSym=0,nIrrep-1
    if (iOper(iSym) == iTest) Invariant = .true.
  end do
  if (Invariant) cycle

  ! Compute total mass and center of mass of the molecule, the center
  ! of the RF cavity is at origin with infinite mass. Hence, the latter
  ! is ignored!

  COM_xyz = Zero
  do iCent=1,nCent
    COM_xyz = COM_xyz+currXYZ(ixyz,iCent)*xMass(iCent)
  end do
  COM_xyz = COM_xyz/TMass

  if (.not. VarT) cycle

  iDeg = 1
  Deg = sqrt(real(iDeg,kind=wp))

  nq = nq+1
  if (.not. Process) mB_Tot = mB_Tot+mB
  if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2
  nqRF = nqRF+1
  write(LuIC,'(A,I2.2,2A)') 'TR',nqRF,' = ',TR_type(ixyz)
  Label = ' '
  write(Label,'(A,I2.2)') 'TR',nqRF

  Val = COM_xyz

  ! Compute the gradient

  Grad(:,:) = Zero
  do iCent=1,nCent
    iAtom = Ind(iCent)
    !write(u6,*) 'iAtom,iCOM=',iAtom,iCOM
    Grad(ixyz,iCent) = dMass(iAtom)/TMass
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('Grad (Trans)',' ',Grad,3,nCent)
# endif

  ! Second derivative is trivially zero!

  Hess(:,:) = Zero
  if (Process) then

    Indq(1,nq) = -2**(ixyz-1)
    Indq(2,nq) = 0
    Indq(3,nq) = 0

    !fconst(nq) = sqrt(Fact*Trans_Const)
    fconst(nq) = sqrt(Trans_Const)
    rMult(nq) = Deg

    Valu(nq,iIter) = Val
    qLbl(nq) = Label

    ! Project the gradient vector

    call ProjSym(nCent,Ind,currXYZ,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq,rMult(nq))

  end if

end do
!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) 'VarR=',VarR
if (VarR) then

  ! A la Malmqvist

  nOrder = 2
  nMass = nCent
  Trans(:) = Zero
  RotVec(:) = Zero
  call mma_allocate(d2RV,3,3*nCent,3*nCent,label='d2RV')
# ifdef _DEBUGPRINT_
  call RecPrt('xMass',' ',xMass,1,nMass)
# endif
  call RotDer(nMass,xMass,currXYZ,ref123,trans,RotAng,RotVec,RotMat,nOrder,dRVdXYZ,d2RV)
# ifdef _DEBUGPRINT_
  call RecPrt('RotVec',' ',RotVec,1,3)
  call RecPrt('RotMat',' ',RotMat,3,3)
  call RecPrt('dRVdXYZ',' ',dRVdXYZ,3,3*nMass)
# endif

  do ixyz=1,3

    Invariant = .false.
    if (ixyz == 1) then
      iTest = 6
    else if (ixyz == 2) then
      iTest = 5
    else
      iTest = 3
    end if
    do iSym=0,nIrrep-1
      if (iOper(iSym) == iTest) Invariant = .true.
    end do
    if (Invariant) cycle

    jxyz = ixyz+1
    if (jxyz > 3) jxyz = 1
    kxyz = jxyz+1
    if (kxyz > 3) kxyz = 1
    iDeg = 1
    Deg = sqrt(real(iDeg,kind=wp))

    nq = nq+1
    if (.not. Process) mB_Tot = mB_Tot+mB
    if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2
    nqRF = nqRF+1
    write(LuIC,'(A,I2.2,2A)') 'TR',nqRF,' = ',TR_type(ixyz+3)
    Label = ' '
    write(Label,'(A,I2.2)') 'TR',nqRF

    Val = RotVec(ixyz)

    ! Compute the gradient

    Grad(:,:) = dRVdXYZ(ixyz,:,:)
#   ifdef _DEBUGPRINT_
    call RecPrt('Grad (Rot)',' ',Grad,3,nCent)
#   endif

    ! Second derivative

    Hess(:,:) = d2RV(ixyz,:,:)

    if (Process) then

      Indq(1,nq) = -(2**(jxyz-1)+2**(kxyz-1))
      Indq(2,nq) = 0
      Indq(3,nq) = 0

      fconst(nq) = sqrt(Rot_Const)
      rMult(nq) = Deg

      Valu(nq,iIter) = Val
      qLbl(nq) = Label

      ! Project the gradient vector

      call ProjSym(nCent,Ind,currXYZ,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq,rMult(nq))

    end if

  end do
  call mma_deallocate(d2RV)
  !write(u6,*) 'nqRF=',nqRF
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(currXYZ)
call mma_deallocate(Ref123)
call mma_deallocate(Grad)
call mma_deallocate(dRVdxyz)
call mma_deallocate(xMass)
call mma_deallocate(Ind)
call mma_deallocate(iDCR)
call mma_deallocate(Hess)

return

end subroutine RF_Coord
