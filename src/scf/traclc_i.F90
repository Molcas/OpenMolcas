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
!               1996, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine TraClc_i(iterLw,nD)
!***********************************************************************
!                                                                      *
! purpose: compute traces                                              *
!                                                                      *
! input:                                                               *
!   Dens    : a few last density matrix differences (nDT,NumDT)        *
!   TwoHam  : a few last two-electron hamiltonians (nDT,NumDT)         *
!   OneHam  : one-electron hamiltonian of length nDT                   *
!   iterLw  : lowest iteration count ...                               *
!             traces are computed from iterLw...iter                   *
!                                                                      *
! output:                                                              *
!   TrDh    : Traces of dD(i)*h of size (nTr)                          *
!   TrDP    : Traces of dD(i)*dP(j) of size (nTr,nTr)                  *
!   TrDD    : Traces of dD(i)*dD(j) of size (nTr,nTr)                  *
!                                                                      *
! called from: OptClc                                                  *
!                                                                      *
! calls to: RWDTG                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Written by:                                                          *
! P.O. Widmark, M.P. Fuelscher and P. Borowski                         *
! University of Lund, Sweden, 1992                                     *
! modified by M. Schuetz, 1996                                         *
! - traces are recomputed for all iterations between iterLw & iter     *
!                                                                      *
!***********************************************************************

use InfSCF, only: Dens, iDisk, iDKeep, Iter, MapDns, nBT, OneHam, TrDD, TrDh, TrDP, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IterLw, nD
integer(kind=iwp) :: i, iD, ii, iPos, iPosL
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iC, iR
#endif
real(kind=wp), allocatable, target :: Aux1(:,:), Aux2(:,:), Aux3(:,:)
real(kind=wp), pointer :: pDens(:,:), pTwoHam(:,:), pVxc(:,:)
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

if (iDKeep < 0) return

!----------------------------------------------------------------------*
!                                                                      *
! Expand the DFT contribution around the reference density             *
! i.e. the last density.                                               *
!                                                                      *
! Note that for ii larger than 1 that the stored density is the        *
! density difference between the two last iterations.                  *
!                                                                      *
!----------------------------------------------------------------------*

! Loop over densities to interpolate over

do ii=iterLw,iter
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Get the one-particle density_ii
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Get pointer to the density

  iPosL = MapDns(ii)

  if (iPosL <= 0) then

    ! If not in memory pick up from disk

    call mma_allocate(Aux1,nBT,nD,Label='Aux1',safe='*')

    ! Pick up the density matrix and external potential

    call RWDTG(-iPosL,Aux1,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    pDens => Aux1
  else
    pDens => Dens(1:nBT,1:nD,iPosL)
  end if

  do iD=1,nD

    ! Trace the one-electron density with the one-electron Hamiltonian.

    TrDh(ii,ii,iD) = DDot_(nBT,pDens(:,iD),1,OneHam,1)

  end do ! iD

  nullify(pDens)
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
end do ! ii
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *

#ifdef _DEBUGPRINT_
do iD=1,nD
  write(u6,'(a)') 'traclc: TrDh'
  write(u6,'(6f16.8)') (TrDh(ii,ii,iD),ii=1,iter)
end do
#endif
call mma_deallocate(Aux1,safe='*')
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
do ii=iterLw,iter

  iPosL = MapDns(ii)
  do iD=1,nD

    ! Trace the two-electron contribution of the Fock matrix with
    ! the density. Diagonal terms.

    if (iPosL > 0) then
      TrDP(ii,ii,iD) = DDot_(nBT,Dens(1,iD,iPosL),1,TwoHam(1,iD,iPosL),1)+DDot_(nBT,Dens(1,iD,iPosL),1,Vxc(1,iD,iPosL),1)
      TrDD(ii,ii,iD) = DDot_(nBT,Dens(1,iD,iPosL),1,Dens(1,iD,iPosL),1)
    else
      write(u6,'(a)') 'traclc: should not happen!!!'
      TrDP(ii,ii,iD) = Zero
      call Abend()
    end if
  end do ! iD

  do i=1,ii-1

    ! Get pointer to density

    iPos = MapDns(i)
    if (iPos <= 0) then
      if (.not. allocated(Aux1)) then
        call mma_allocate(Aux1,nBT,nD,Label='Aux1')
        call mma_allocate(Aux2,nBT,nD,Label='Aux2')
        call mma_allocate(Aux3,nBT,nD,Label='Aux3')
      end if
      call RWDTG(-iPos,Aux1,nBT*nD,'R','TWOHAM',iDisk,size(iDisk,1))
      call RWDTG(-iPos,Aux2,nBT*nD,'R','dVxcdR',iDisk,size(iDisk,1))
      call RWDTG(-iPos,Aux3,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
      pTwoHam => Aux1
      pVxc => Aux2
      pDens => Aux3
    else
      pTwoHam => TwoHam(1:nBT,1:nD,iPos)
      pVxc => Vxc(1:nBT,1:nD,iPos)
      pDens => Dens(1:nBT,1:nD,iPos)
    end if

    do iD=1,nD
      TrDP(i,ii,iD) = DDot_(nBT,Dens(1,iD,iPosL),1,pTwoHam(:,iD),1)
      TrDP(ii,i,iD) = TrDP(i,ii,iD)
      TrDP(i,ii,iD) = TrDP(i,ii,iD)+DDot_(nBT,Dens(1,iD,iPosL),1,pVxc(:,iD),1)
      TrDP(ii,i,iD) = TrDP(ii,i,iD)+DDot_(nBT,Vxc(1,iD,iPosL),1,pDens(:,iD),1)
      TrDD(i,ii,iD) = DDot_(nBT,Dens(1,iD,iPosl),1,pDens(:,iD),1)
      TrDD(ii,i,iD) = TrDD(i,ii,iD)

    end do ! iD

#   ifdef _DEBUGPRINT_
    do iD=1,nD
      write(u6,*) 'iteration:',ii
      write(u6,'(a)') 'traclc: TrDh'
      do iR=1,ii
        write(u6,'(6f16.8)') TrDh(iR,iR,iD)
      end do
      write(u6,'(a)') 'traclc: TrDP'
      do iR=1,ii
        write(u6,'(6f16.8)') (TrDP(iR,iC,iD),iC=1,ii)
      end do
      write(u6,'(a)') 'traclc: TrDD'
      do iR=1,ii
        write(u6,'(6f16.8)') (TrDD(iR,iC,iD),iC=1,ii)
      end do
    end do ! iD
#   endif
    nullify(pTwoHam,pVxc,pDens)

  end do ! i

end do ! ii

if (allocated(Aux1)) then
  call mma_deallocate(Aux1)
  call mma_deallocate(Aux2)
  call mma_deallocate(Aux3)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

#undef _DEBUGPRINT_
end subroutine TraClc_i
