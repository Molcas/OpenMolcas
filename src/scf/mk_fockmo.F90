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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Mk_FockMO(O,S,nOTSD,C,nC,FockMO,nFockMO,nD,iOpt)

use Index_Functions, only: nTri_Elem
use InfSCF, only: iDisk, MapDns, MaxBas, nBas, nBO, nBT, nFro, nnFr, nOcc, nOrb, nSym, OrbType, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nOTSD, nC, nFockMO, nD, iOpt
real(kind=wp), intent(in) :: O(nOTSD), S(nOTSD), C(nC,nD)
real(kind=wp), intent(out) :: FockMO(nFockMO,nD)
integer(kind=iwp) :: i, iD, ig, ih, ij, iOff, iSym, it, j, jDT, k, l, nBs, nOr, nOrbmF
real(kind=wp), allocatable :: Aux1(:), Aux2(:), FckM(:,:)
real(kind=wp), allocatable, target :: AuxT(:,:), AuxV(:,:)
real(kind=wp), pointer :: T(:,:), V(:,:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
jDT = MapDns(iOpt)
if (jDT < 0) then
  call mma_allocate(AuxT,nOTSD,nD,Label='AuxT')
  call mma_allocate(AuxV,nOTSD,nD,Label='AuxV')

  call RWDTG(-jDT,AuxT,nOTSD*nD,'R','TWOHAM',iDisk,size(iDisk,1))
  call RWDTG(-jDT,AuxV,nOTSD*nD,'R','dVxcdR',iDisk,size(iDisk,1))

  T(1:nOTSD,1:nD) => AuxT(:,:)
  V(1:nOTSD,1:nD) => AuxV(:,:)
else
  T(1:nOTSD,1:nD) => TwoHam(:,:,jDT)
  V(1:nOTSD,1:nD) => Vxc(:,:,jDT)
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'Mk_FockMO: input arrays'
write(u6,*) '==================================================='
call NrmClc(O,nOTSD,'EGrad','O')
call NrmClc(T,nOTSD*nD,'EGrad','T')
call NrmClc(V,nOTSD*nD,'EGrad','V')
call NrmClc(C,nC*nD,'EGrad','C')
!do iD=1,nD
!  write(u6,*) 'OrbType(:,iD)',OrbType(:,iD)
!end do ! iD
write(u6,*) '==================================================='
write(u6,*)
#endif
!----------------------------------------------------------------------*

! Allocate memory for modified fock matrix
call mma_allocate(FckM,nBT,nD,Label='FckM')
FckM(:,:) = Zero
FockMO(:,:) = Zero

! Allocate memory for auxiliary matrices
call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
call mma_allocate(Aux2,MaxBas**2,Label='Aux2')

do iD=1,nD

  FckM(:,iD) = O(:)+T(:,iD)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iD=',iD
  call NrmClc(FckM(1,iD),nBT,'Mk_FockMO','FckM')
# endif
  if (nnFr > 0) call ModFck(FckM(:,iD),S,nBT,C(:,iD),nBO,nOcc(:,1))

  FckM(:,iD) = FckM(:,iD)+V(:,iD)
# ifdef _DEBUGPRINT_
  call NrmClc(FckM(1,iD),nBT,'Mk_FockMO','FckM')
# endif

  iOff = 0
  ij = 1
  it = 1
  iG = 1
  do iSym=1,nSym
    nBs = nBas(iSym)
    nOr = nOrb(iSym)
    nOrbmF = nOrb(iSym)-nFro(iSym)

    if (nOrb(iSym) > 0) then

      ! Square Fock matrix and perform C(T)F
      Aux2(:) = Zero
      call Square(FckM(ij,iD),Aux2,1,nBs,nBs)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iSym=',iSym
      call NrmClc(Aux2,nBs*nBs,'Mk_FockMO','Aux2')
#     endif
      Aux1(:) = Zero
      call DGEMM_('T','N',nOr,nBs,nBs, &
                  One,C(it,iD),nBs, &
                  Aux2,nBs, &
                  Zero,Aux1,nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux1,nOr*nBs,'Mk_FockMO','Aux1')
#     endif
      ! C(T)FDSC
      Aux2(:) = Zero
      call DGEMM_('N','N',nOr,nOr,nBs, &
                  One,Aux1,nOr, &
                  C(it,iD),nBs, &
                  Zero,FockMO(iG,iD),nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(FockMO(iG,iD),nOr*nOr,'Mk_FockMO','FockMO')
#     endif

      ! At this point enforce that the gradient is exactly zero
      ! for elements corresponding to orbitals of different
      ! fermion types.

      do i=1,nOr
        if (i <= nFro(iSym)) then
          k = -1
        else
          k = OrbType(iOff+i-nFro(iSym),iD)
        end if

        do j=1,nOr
          if (j <= nFro(iSym)) then
            l = -1
          else
            l = OrbType(iOff+j-nFro(iSym),iD)
          end if

          ih = iG+(i-1)*nOr+j-1
          if (k /= l) FockMO(ih,iD) = Zero

        end do
      end do

    end if
    ij = ij+nTri_Elem(nBs)
    it = it+nBs*nOr
    iG = iG+nOr*nOr
    iOff = iOff+nOrbmF

  end do ! iSym

end do ! iD

! Deallocate memory
call mma_deallocate(Aux2)
call mma_deallocate(Aux1)
call mma_deallocate(FckM)

if (jDT < 0) then
  call mma_deallocate(AuxT)
  call mma_deallocate(AuxV)
end if
nullify(T,V)
#ifdef _DEBUGPRINT_
call RecPrt('Mk_FockMO: FockMO',' ',FockMO,nFockMO,nD)
#endif

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine Mk_FockMO
