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
subroutine Fix_Ref(Old_Ref)
! Transform coordinates, displacements and gradients
! to make them consistent with a new reference

use LnkLst, only: GetNod, iVPtr, LLDelt, LLdGrd, LLGrad, LLlGrd, LLx, LstPtr, PutVec, SCF_V
use InfSCF, only: Iter, Iter_Start, kOV, mOV, nD, nFro, nOcc, nOFS, nOrb, nSym, TimFld
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Old_Ref
integer(kind=iwp) :: i, iD, iEnd, inode, iORef, iO, iOff1, iOff2, iSt, iSyBlpt, iSym, nO, nOccmF(8,nD), nOF, nVrt
real(kind=wp) :: Cpu1, Cpu2, Tim1, Tim2, Tim3
logical(kind=iwp) :: NotUnit
real(kind=wp), allocatable :: dG(:), dX(:), G(:), RedRot(:,:), RoM(:), RotRef(:,:), tmp(:), X(:), XRef(:)
logical(kind=iwp), external :: isUnit
real(kind=wp), parameter :: thrs = 1.0e-10

call Timing(Cpu1,Tim1,Tim2,Tim3)

call mma_allocate(X,mOV,Label='X')
call mma_allocate(dX,mOV,Label='dX')
call mma_allocate(G,mOV,Label='G')
call mma_allocate(dG,mOV,Label='dG')
call mma_allocate(RoM,nOFS,Label='RoM')
call mma_allocate(RedRot,nOFS,nD,Label='RedRot')
call mma_allocate(tmp,mOV,Label='tmp')

iORef = LstPtr(Old_Ref,LLx)   ! Pointer to X_old(i_oldref)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Iter=',Iter
write(u6,*) 'Iter_Start=',Iter_Start
write(u6,*) 'Old_Ref=',Old_Ref
call NrmClc(SCF_V(iORef)%A,mOV,'Fix_Ref','X(i_oldref)')
write(u6,*)
#endif

! Get the displacement for the new reference
call mma_allocate(XRef,mOV,Label='XRef')
XRef(:) = -SCF_V(iORef)%A(:)
! Get the rotation matrix for the new reference
call mma_allocate(RotRef,nOFS,nD,Label='RotRef')
iEnd = 0
do iD=1,nD
  iSt = iEnd+1
  iEnd = iEnd+kOV(iD)
  nOccmF(:,iD) = nOcc(:,iD)-nFro(:)
  call ExpKap(XRef(iSt:iEnd),kOV(iD),RotRef(:,iD),nOccmF(:,iD))
end do

! Loop over all iterations starting at Iter_Start

do i=Iter_Start,Iter

  call GetNod(i,LLx,inode)
  if (inode == 0) then
    write(u6,*) 'inode == 0'
    call Abend()
  end if
  call iVPtr(X,mOV,inode)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'X(i) before  i=',i
  call NrmClc(X,mOV,'Fix_Ref','X')
# endif

  !***************************************************************
  ! Transform the coordinates (kappa matrix, rotation parameters)
  !***************************************************************

  X(:) = X(:)+XRef(:)

  iEnd = 0
  do iD=1,nD
    iSt = iEnd+1
    iEnd = iEnd+kOV(iD)
    iOff1 = iSt-1
    ! Get the rotation matrix from the old reference
    call ExpKap(X(iSt:iEnd),kOV(iD),RoM,nOccmF(:,iD))
    iSyBlpt = 1
    do iSym=1,nSym
      nOF = nOrb(iSym)-nFro(iSym)
      nVrt = nOrb(iSym)-nOcc(iSym,iD)
      nO = nOccmF(iSym,iD)
      if (nVrt*nO < 1) cycle
      ! Calculate the combined rotation matrix from the new reference
      call dgemm_('T','N',nOF,nOF,nOF,One,RotRef(iSyBlpt,iD),nOF,RoM(iSyBlpt),nOF,Zero,RedRot(iSyBlpt,iD),nOF)
      ! Obtain the kappa matrix and redundant rotations corresponding to the combined rotation
      call log_SVD(nOF,nO,RedRot(iSyBlpt,iD))
      iOff2 = iSyBlpt+nO-1
      do iO=1,nO
        X(iOff1+1:iOff1+nVrt) = RedRot(iOff2+1:iOff2+nVrt,iD)
        iOff1 = iOff1+nVrt
        iOff2 = iOff2+nOF
      end do
      iSyBlpt = iSyBlpt+nOF**2
    end do
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'X(i) after  i=',i
  call NrmClc(X,mOV,'Fix_Ref','X')
# endif

  call PutVec(X,mOV,i,'OVWR',LLx)

  !*****************************
  ! Recompute the displacements
  ! dX(i-1) = X(i)-X(i-1)
  !*****************************

  if (i > Iter_Start) then
    dX(:) = X(:)-dX(:)
    call PutVec(dX,mOV,i-1,'OVWR',LLDelt)
  end if
  dX(:) = X(:)

  !*****************************************************
  ! Transform the gradients
  ! (these are local gradients at the rotated orbitals)
  !*****************************************************

  call GetNod(i,LLlGrd,inode)
  if (inode == 0) then
    write(u6,*) 'inode == 0'
    call Abend()
  end if
  call iVPtr(G,mOV,inode)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'G(i) before  i=',i
  call NrmClc(G,mOV,'Fix_Ref','G')
# endif

  iEnd = 0
  do iD=1,nD
    iSt = iEnd+1
    iEnd = iEnd+kOV(iD)
    iOff1 = iSt
    iSyBlpt = 1
    do iSym=1,nSym
      nOF = nOrb(iSym)-nFro(iSym)
      nVrt = nOrb(iSym)-nOcc(iSym,iD)
      nO = nOccmF(iSym,iD)
      if (nVrt*nO < 1) cycle
      ! This is a straightforward multiplication by the redundant rotations
      ! Only transform if not a unit matrix
      iOff2 = iSyBlpt+nO*(nOF+1)
      if (isUnit(RedRot(iSyBlpt,iD),nO,nOF,thrs)) then
        NotUnit = .not. isUnit(RedRot(iOff2,iD),nVrt,nOF,thrs)
      else
        NotUnit = .true.
      end if
      if (NotUnit) then
        call dgemm_('N','T',nVrt,nO,nO,One,G(iOff1),nVrt,RedRot(iSyBlpt,iD),nOF,Zero,tmp,nVrt)
        call dgemm_('N','N',nVrt,nO,nVrt,One,RedRot(iOff2,iD),nOF,tmp,nVrt,Zero,G(iOff1),nVrt)
      end if
      iOff1 = iOff1+nO*nVrt
      iSyBlpt = iSyBlpt+nOF**2
    end do
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'G(i) after  i=',i
  call NrmClc(G,mOV,'Fix_Ref','loc G')
# endif

  call PutVec(G,mOV,i,'OVWR',LLlGrd)

  !*****************************************************
  ! Transform the gradients to the new reference
  !*****************************************************

  iEnd = 0
  do iD=1,nD
    iSt = iEnd+1
    iEnd = iEnd+kOV(iD)
    call TrGrad(X(iSt),kOV(iD),G(iSt),nOccmF(:,iD))
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'G(i) after  i=',i
  call NrmClc(G,mOV,'Fix_Ref','G')
# endif

  call PutVec(G,mOV,i,'OVWR',LLGrad)

  !************************************
  ! Recompute the gradient differences
  ! dG(i-1) = G(i)-G(i-1)
  !************************************

  if (i > Iter_Start) then
    dG(:) = G(:)-dG(:)
    call PutVec(dG,mOV,i-1,'OVWR',LLdGrd)
  end if
  dG(:) = G(:)

end do

call mma_deallocate(RotRef)

call mma_deallocate(X)
call mma_deallocate(dX)
call mma_deallocate(G)
call mma_deallocate(dG)
call mma_deallocate(RoM)
call mma_deallocate(RedRot)
call mma_deallocate(tmp)
call mma_deallocate(XRef)

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(11) = TimFld(11)+(Cpu2-Cpu1)

return

end subroutine Fix_Ref
