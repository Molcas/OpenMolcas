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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2017,2022, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine EGrad(O,S,nOTSD,C,nC,G,nG,nD,iOpt)
!***********************************************************************
!                                                                      *
!     purpose: This routine calculates the gradient of the SCF energy  *
!              with respect to the rotation parameters.                *
!              Grad(E) = C(t)(FDS-SDF)C                                *
!                                                                      *
!                                                                      *
!     input:                                                           *
!       O       : one-electron hamiltonian of length nOTSD             *
!       S       : overlap in AO basis of length nOTSD                  *
!       C       : matrix transforming to the set of orthonormal        *
!                 (and spherical, if needed) functions of length nC    *
!                                                                      *
!     output:                                                          *
!       G       : gradient of the SCF energy with respect to the       *
!                 rotation parameters of length nG                     *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: Dens, iDisk, MapDns, MaxBas, nBas, nBO, nBT, nFro, nnFr, nOcc, nOrb, nSym, OrbType, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nOTSD, nC, nG, nD, iOpt
real(kind=wp), intent(in) :: O(nOTSD), S(nOTSD), C(nC,nD)
real(kind=wp), intent(inout) :: G(nG,nD)
integer(kind=iwp) :: i, iD, ig, ih, ij, iOff, iSym, it, j, jDT, k, l, nBs, nOr, nOrbmF
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i_Max, j_Max
real(kind=wp) :: GMax
#endif
real(kind=wp), allocatable :: Aux1(:), Aux2(:), Aux3(:), FckM(:,:)
real(kind=wp), allocatable, target :: AuxD(:,:), AuxT(:,:), AuxV(:,:)
real(kind=wp), pointer :: D(:,:), T(:,:), V(:,:)

!----------------------------------------------------------------------*
! Start
jDT = MapDns(iOpt)
if (jDT < 0) then
  call mma_allocate(AuxD,nOTSD,nD,Label='AuxD')
  call mma_allocate(AuxT,nOTSD,nD,Label='AuxT')
  call mma_allocate(AuxV,nOTSD,nD,Label='AuxV')

  call RWDTG(-jDT,AuxD,nOTSD*nD,'R','DENS  ',iDisk,size(iDisk,1))
  call RWDTG(-jDT,AuxT,nOTSD*nD,'R','TWOHAM',iDisk,size(iDisk,1))
  call RWDTG(-jDT,AuxV,nOTSD*nD,'R','dVxcdR',iDisk,size(iDisk,1))

  D(1:nOTSD,1:nD) => AuxD(:,:)
  T(1:nOTSD,1:nD) => AuxT(:,:)
  V(1:nOTSD,1:nD) => AuxV(:,:)
else
  D(1:nOTSD,1:nD) => Dens(:,:,jDT)
  T(1:nOTSD,1:nD) => TwoHam(:,:,jDT)
  V(1:nOTSD,1:nD) => Vxc(:,:,jDT)
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'EGrad: input arrays'
write(u6,*) '==================================================='
call NrmClc(O,nOTSD,'EGrad','O')
call NrmClc(S,nOTSD,'EGrad','S')
call NrmClc(D,nOTSD*nD,'EGrad','D')
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
G(:,:) = Zero

! Allocate memory for auxiliary matrices
call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
call mma_allocate(Aux2,MaxBas**2,Label='Aux2')
call mma_allocate(Aux3,MaxBas**2,Label='Aux3')

do iD=1,nD

  FckM(:,iD) = O(:)+T(:,iD)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iD=',iD
  call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
# endif
  if (nnFr > 0) call ModFck(FckM(:,iD),S,nBT,C(:,iD),nBO,nOcc(:,iD))

  FckM(:,iD) = FckM(:,iD)+V(:,iD)
# ifdef _DEBUGPRINT_
  call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
# endif

  iOff = 0
  ij = 1
  it = 1
  ig = 1
  do iSym=1,nSym
    nBs = nBas(iSym)
    nOr = nOrb(iSym)
    nOrbmF = nOrb(iSym)-nFro(iSym)

    if (nOrb(iSym) > 0) then

      ! Square Fock matrix and perform C(T)F
      Aux2(:) = Zero
      call Square(FckM(ij:,iD),Aux2,1,nBs,nBs)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iSym=',iSym
      call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#     endif
      Aux1(:) = Zero
      call DGEMM_('T','N',nOr,nBs,nBs, &
                  One,C(it,iD),nBs, &
                  Aux2,nBs, &
                  Zero,Aux1,nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux1,nOr*nBs,'EGrad','Aux1')
#     endif

      ! Square density matrix and perform C(T)FD
      Aux2(:) = Zero
      call DSq(D(ij:,iD),Aux2,1,nBs,nBs)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#     endif
      Aux3(:) = Zero
      call DGEMM_('N','N',nOr,nBs,nBs, &
                  One,Aux1,nOr, &
                  Aux2,nBs, &
                  Zero,Aux3,nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux3,nOr*nBs,'EGrad','Aux3')
#     endif

      ! Square overlap matrix and perform C(T)FDS
      Aux2(:) = Zero
      call Square(S(ij:),Aux2,1,nBs,nBs)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#     endif
      Aux1(:) = Zero
      call DGEMM_('N','N',nOr,nBs,nBs, &
                  One,Aux3,nOr, &
                  Aux2,nBs, &
                  Zero,Aux1,nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux1,nOr*nBs,'EGrad','Aux1')
#     endif
      ! C(T)FDSC
      Aux2(:) = Zero
      call DGEMM_('N','N',nOr,nOr,nBs, &
                  One,Aux1,nOr, &
                  C(it,iD),nBs, &
                  Zero,Aux2,nOr)
#     ifdef _DEBUGPRINT_
      call NrmClc(Aux2,nOr*nOr,'EGrad','Aux2')
#     endif

      call Asym(Aux2,G(ig:ig+nOr**2,iD),nOr)
#     ifdef _DEBUGPRINT_
      write(u6,*)
      call NrmClc(G,nG*nD,'EGrad','G')
      write(u6,*)
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

          ih = ig+(i-1)*nOr+j-1
          if (k /= l) G(ih,iD) = Zero

        end do
      end do

    end if
    ij = ij+nTri_Elem(nBs)
    it = it+nBs*nOr
    ig = ig+nOr**2
    iOff = iOff+nOrbmF

  end do ! iSym

end do ! iD

! Deallocate memory
call mma_deallocate(Aux3)
call mma_deallocate(Aux2)
call mma_deallocate(Aux1)
call mma_deallocate(FckM)

G(:,:) = Two*G(:,:)

#ifdef _DEBUGPRINT_
i_Max = 0
j_Max = 0
GMax = Zero
do i=1,nD
  do j=1,nG
    if (GMax < abs(G(j,i))) then
      GMax = abs(G(j,i))
      i_Max = i
      j_Max = j
    end if
  end do
end do
write(u6,*) 'GMax,i_Max,j_Max=',GMax,i_Max,j_Max
call NrmClc(G,nG*nD,'EGrad','G')
#endif
if (jDT < 0) then
  call mma_deallocate(AuxD)
  call mma_deallocate(AuxT)
  call mma_deallocate(AuxV)
end if
nullify(D,T,V)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

contains

subroutine Asym(H,A,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: H(n,n)
  real(kind=wp), intent(inout) :: A(n,n)
  integer(kind=iwp) :: i

  do i=1,n
    A(i,1:i-1) = H(i,1:i-1)-H(1:i-1,i)
    A(i,i) = Zero
  end do

end subroutine ASym

end subroutine EGrad
