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

subroutine CI_KAP(ipcid,fock,fockOut,isym)

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: ipnout !, W
use MCLR_Data, only: ipCI, n2Dens, nDens, nNA
use input_mclr, only: ntAsh, State_Sym
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipCID, iSym
real(kind=wp), intent(out) :: Fock(nDens), FockOut(nDens)
integer(kind=iwp) :: i, ij, ij1, j, k, kl1, l, nDim
real(kind=wp) :: D0
real(kind=wp), allocatable :: De(:), Pe(:), tmpDe(:,:), tmpDeM(:,:), tmpP(:), tmpPM(:,:,:,:)

call ipnout(-1)
call mma_allocate(De,ntash**2,Label='De')
call mma_allocate(Pe,n2Dens,Label='Pe')

call CIDens_SA(.true.,ipCI,ipCid,State_sym,State_Sym,Pe,De)

d0 = Zero
! ======================================================================
if (doDMRG) then
  call dmrg_dim_change_mclr(LRras2,ntash,0)
  call dmrg_dim_change_mclr(RGras2,ndim,0)

  call mma_allocate(tmpDe,ndim,ndim,Label='tmpDe')
  call mma_allocate(tmpP,nTri_Elem(ndim**2),Label='tmpP')
  call mma_allocate(tmpDeM,ntash,ntash,Label='tmpDeM')
  call mma_allocate(tmpPM,ntash,ntash,ntash,ntash,Label='tmpPM')
  tmpDe(:,:) = Zero
  tmpP(:) = Zero
  tmpDeM(:,:) = Zero
  tmpPM(:,:,:,:) = Zero

  ij = 0
  do i=1,ntash
    do j=1,ntash
      ij = ij+1
      if (abs(De(ij)) < 1.0e-12_wp) De(ij) = Zero
      tmpDeM(i,j) = De(ij)
    end do
  end do

  ij = 0
  do i=1,ndim
    do j=1,ndim
      ij = ij+1
      if ((i > ntash) .or. (j > ntash)) then
        tmpDe(i,j) = Zero
      else
        tmpDe(i,j) = tmpDeM(i,j)
      end if
    end do
  end do

  do i=1,ntash
    do j=1,ntash
      do k=1,ntash
        do l=1,ntash
          ij1 = ntash*(i-1)+j
          kl1 = ntash*(k-1)+l
          if (ij1 >= kl1) then
            if (abs(Pe(iTri(ij1,kl1))) < 1.0e-12_wp) Pe(iTri(ij1,kl1)) = Zero
            tmpPM(i,j,k,l) = Pe(iTri(ij1,kl1))
          end if
        end do
      end do
    end do
  end do

  do i=1,ndim
    do j=1,ndim
      do k=1,ndim
        do l=1,ndim
          ij1 = ndim*(i-1)+j
          kl1 = ndim*(k-1)+l
          if (ij1 >= kl1) then
            if ((i > ntash) .or. (j > ntash) .or. (k > ntash) .or. (l > ntash)) then
              tmpP(iTri(ij1,kl1)) = Zero
            else
              tmpP(iTri(ij1,kl1)) = tmpPM(i,j,k,l)
            end if
          end if
        end do
      end do
    end do
  end do

  ! ====================================================================
  call dmrg_dim_change_mclr(RGras2,ntash,0)
  call dmrg_dim_change_mclr(RGras2,nna,0)

  !call ipin(ipCID)
  !call ipin(ipci)
  !call projecter(W(ipCID)%A,W(ipci)%A,De,Pe)
  call FockGen(d0,tmpDe,tmpP,Fock,FockOut,isym) ! yma modified

else

  !call ipin(ipCID)
  !call ipin(ipci)
  !call projecter(W(ipCID)%A,W(ipci)%A,De,Pe)
  call FockGen(d0,De,Pe,Fock,FockOut,isym)
end if
call mma_deallocate(De)
call mma_deallocate(Pe)

! ======================================================================
if (doDMRG) then  !yma
  call dmrg_dim_change_mclr(LRras2,nna,0)
  call dmrg_dim_change_mclr(LRras2,ntash,0)
  call mma_deallocate(tmpDe)
  call mma_deallocate(tmpDeM)
  call mma_deallocate(tmpP)
  call mma_deallocate(tmpPM)
end if
! ======================================================================

end subroutine CI_KAP
