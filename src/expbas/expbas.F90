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
! Copyright (C) 2008, Bjorn O. Roos                                    *
!               2008, Valera Veryazov                                  *
!***********************************************************************

subroutine expbas(ireturn)
!***********************************************************************
!                                                                      *
!     Objective: Expand MOs to larger basis set                        *
!                                                                      *
!     B. O. Roos, University of Lund, April 2008.                      *
!                                                                      *
!***********************************************************************

use info_expbas_mod, only: EB_FileOrb, LenIn, mxsym, n_orb_kinds
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, ib1, ib2, iErr, iLen, ind, ishift, ist1, ist2, iSym, Lu_, LuInpOrb, nb1, nb2, nDim1, nDim2, nSym1, nSym2, &
                     nTot1, nTot2, nBas1(mxsym), nBas2(mxsym)
character(len=80) :: VecTit
character(len=512) :: FName
logical(kind=iwp) :: Exist_1, Exist_2, okay
integer(kind=iwp), allocatable :: indt1(:), indt2(:), IndType(:)
real(kind=wp), allocatable :: CMO1(:), CMO2(:), Eorb1(:), Eorb2(:), Occ1(:), Occ2(:)
character(len=LenIn+8), allocatable :: Bas1(:), Bas2(:)

!----------------------------------------------------------------------*
!     Read information from Runfile 1                                  *
!----------------------------------------------------------------------*
FName = 'RUNFIL1'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),Exist_1)
if (.not. Exist_1) then
  write(u6,*) 'Error finding file '//FName(:iLen)
  call Abend()
end if
call namerun(FName(:iLen))
call get_iScalar('nSym',nSym1)
call Get_iArray('nBas',nBas1,nSym1)
nDim1 = 0
nTot1 = 0
do iSym=1,nSym1
  nDim1 = nDim1+nBas1(iSym)
  nTot1 = nTot1+nBas1(iSym)**2
end do
call mma_allocate(Bas1,nDim1,label='Bas1')
call Get_cArray('Unique Basis Names',Bas1,(Lenin+8)*nDim1)
!----------------------------------------------------------------------*
!     Read information from Runfile 2                                  *
!----------------------------------------------------------------------*
FName = 'RUNFIL2'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),Exist_2)
if (.not. Exist_2) then
  write(u6,*) 'Error finding file '//FName(:iLen)
  call Abend()
end if
call namerun(FName(:iLen))
call get_iScalar('nSym',nSym2)
call Get_iArray('nBas',nBas2,nSym2)
nDim2 = 0
nTot2 = 0
do iSym=1,nSym2
  nDim2 = nDim2+nBas2(iSym)
  nTot2 = nTot2+nBas2(iSym)**2
end do
call mma_allocate(Bas2,nDim2,label='Bas2')
call Get_cArray('Unique Basis Names',Bas2,(LenIn+8)*nDim2)
!----------------------------------------------------------------------*
!     Read MO coefficients from a formatted vector file                *
!----------------------------------------------------------------------*
call mma_allocate(CMO1,nTot1,label='CMO1')
call mma_allocate(Eorb1,nDim1,label='Eorb1')
call mma_allocate(Occ1,nDim1,label='Occ1')
call mma_allocate(indt1,nDim1,label='indt1')
call mma_allocate(CMO2,nTot2,label='CMO2')
call mma_allocate(Eorb2,nDim2,label='Eorb2')
call mma_allocate(Occ2,nDim2,label='Occ2')
call mma_allocate(indt2,nDim2,label='indt2')
FName = EB_FileOrb
if (len_trim(FName) == 0) FName = 'INPORB'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),okay)
if (okay) then
  LuInpOrb = 50
  call RdVec(FName(:iLen),LuInpOrb,'COEI',nSym1,nBas1,nBas1,CMO1,Occ1,Eorb1,indt1,VecTit,1,iErr)
else
  write(u6,*) 'RdCMO: Error finding MO file'
  call Abend()
end if

!----------------------------------------------------------------------*
!     Print and check input information                                *
!----------------------------------------------------------------------*
!write(u6,910) 'Start of option expand.'
write(u6,930) trim(Vectit)
write(u6,910) 'Information from input runfile'
write(u6,920) 'Number of symmetries',nSym1
write(u6,920) 'Number of basis functions',(nBas1(i),i=1,nSym1)
write(u6,910) 'Information from expanded basis set runfile'
write(u6,920) 'Number of symmetries',nSym2
write(u6,920) 'Number of basis functions',(nBas2(i),i=1,nSym2)
! Check for inconsistensies:
if (nSym1 /= nSym2) then
  write(u6,*) 'Symmetries are not equal. Stop here',nSym1,nSym2
  call Abend()
end if
do isym=1,nSym1
  if (nBas1(isym) > nBas2(isym)) then
    write(u6,*) 'Second basis set must be larger than first'
    write(u6,*) 'not fulfilled in sym',isym,'basis functions are',nBas1(isym),nBas2(isym)
    call Abend()
  end if
end do
!----------------------------------------------------------------------*
!     Build the new orbitals                                           *
!----------------------------------------------------------------------*
ist1 = 1
ist2 = 1
ib1 = 1
ib2 = 1
do isym=1,nsym1
  nb1 = nBas1(isym)
  nb2 = nBas2(isym)
  if (nb2 > 0) then
    call expandbas(Bas1(ib1),nb1,Bas2(ib2),nb2,CMO1(ist1),CMO2(ist2),Occ1(ib1),Eorb1(ib1),indt1(ib1), &
                                                                     Occ2(ib2),Eorb2(ib2),indt2(ib2))
    ist1 = ist1+nb1**2
    ist2 = ist2+nb2**2
    ib1 = ib1+nb1
    ib2 = ib2+nb2
  end if
end do
call mma_deallocate(Bas1)
call mma_deallocate(Bas2)
!----------------------------------------------------------------------*
!     Write the new orbitals in to the file EXPORB                     *
!----------------------------------------------------------------------*
! First resort indt to standard
call mma_allocate(IndType,nSym2*n_orb_kinds,label='IndType')
IndType(:) = 0
ind = 0
ishift = 0
do isym=1,nSym2
  nb2 = nBas2(isym)
  if (nb2 /= 0) then
    do ib2=1,nb2
      ind = ind+1
      IndType(ishift+indt2(ind)) = IndType(ishift+indt2(ind))+1
    end do
  end if
  ishift = ishift+n_orb_kinds
end do

VecTit = 'Basis set expanded orbital file EXPORB'
Lu_ = 60
call WRVEC('EXPORB',LU_,'COEI',nSym2,nBas2,nBas2,CMO2,Occ2,Eorb2,IndType,VecTit)
write(u6,*) 'New orbitals have been built in file EXPORB'
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
call mma_deallocate(CMO1)
call mma_deallocate(Eorb1)
call mma_deallocate(Occ1)
call mma_deallocate(indt1)
call mma_deallocate(CMO2)
call mma_deallocate(Eorb2)
call mma_deallocate(Occ2)
call mma_deallocate(indt2)
call mma_deallocate(IndType)

ireturn = 0

return

910 format(/1x,a)
920 format(1x,a30,8i5)
930 format(/1x,'Header on input orbitals file:'/a)

end subroutine expbas
