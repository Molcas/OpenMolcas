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
! Copyright (C) 2007, Bjorn O. Roos                                    *
!***********************************************************************

subroutine AFREEZ(NSYM,NBAS,NFRO,NISH,NASH,NSSH,NDEL,NAME,nName,NAMFRO,LNFRO,DPQ,nDPQ,THRFR,THRDE,IFQCAN,CMO,NCMO)
!****************************************************************************
!                                                                           *
! Purpose: to select orbitals, which will be frozen in the CASPT2           *
! calculations based on a selection of atoms controlled by the input        *
! keyword AFREeze.                                                          *
! Each inactive orbital is checked for the fraction of electrons located    *
! on the selected atoms. If smaller than a given threshold, the orbital     *
! will be frozen.                                                           *
! Called by READIN_CASPT2                                                   *
! Author: B. O. Roos in July 2007 for MOLCAS-7                              *
!     Calling parameters:                                                   *
!     NSYM   : Number of symmetries                                         *
!     NFRO   : Number of frozen orbitals (modified by the program)          *
!     NISH   : Number of inactive orbitals                                  *
!     Name   : Center and function type label per basis function            *
!     Namfro : names of atoms to be selected (length lnfro)                 *
!     Labfro : labels for orbitals to be frozen                             *
!     CMO    : Orbital coefficients                                         *
!     DPQ    : The charge matrix for a given orbital                        *
!     THRFR : Threshold for freezing orbitals                               *
!     THRDE : Threshold for deleting orbitals                               *
!                                                                           *
!****************************************************************************

use OneDat, only: sNoNuc, sNoOri
use definitions, only: iwp, wp, u6
use Molcas, only: LenIn, MxBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

implicit none
integer(kind=iwp), intent(in) :: NSYM
integer(kind=iwp), intent(in) :: NBAS(NSYM), NASH(NSYM)
integer(kind=iwp), intent(inout) :: NFRO(NSYM), NISH(NSYM), NSSH(NSYM), NDEL(NSYM)
integer(kind=iwp), intent(in) :: nName
character(len=LenIn+8), intent(in) :: NAME(nName)
integer(kind=iwp), intent(in) :: LnFro
character(len=4), intent(in) :: NAMFRO(LnFro)
integer(kind=iwp), intent(in) :: nDPQ
real(kind=wp), intent(out) :: DPQ(nDPQ)
real(kind=wp), intent(in) :: THRFR, THRDE
integer(kind=iwp), intent(inout) :: IFQCAN
integer(kind=iwp), intent(in) :: NCMO
real(kind=wp), intent(inout) :: CMO(nCMO)
integer(kind=iwp) :: LABFRO(mxbas)
real(kind=wp), allocatable :: SMAT(:)
character(len=8) :: Label
real(kind=wp), parameter :: Thrs = 1.0e-6_wp
real(kind=wp) chksum, selch, Swap
integer(kind=iwp) :: I, ib, iComp, imo, imo0, iname, iopt, ipp, ipq, ipq0, iqq, irc, ist1, ist2, isym, isymlbl, nb2, NBAST, nbi, &
                     ndi, nfi, nfro1, ni, np, nq, nsi, NSMAT, ntri, ipq1, nin

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

!write(u6,*) 'Entering AFreez'
NBAST = 0
ntri = 0
do I=1,NSYM
  NBAST = NBAST+NBAS(I)
  ntri = (nbas(i)+nbas(i)**2)/2+ntri
end do
if (NBAST > MXBAS) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
NSMAT = NTRI+6
call MMA_ALLOCATE(SMAT,NSMAT)
isymlbl = 1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'Mltpl  0'
iComp = 1
call RdOne(irc,iopt,Label,iComp,SMAT,isymlbl)

!----------------------------------------------------------------------*
!write(u6,*) 'molecular orbitals before localization'
!imo = 0
!do isym=1,nsym
!  nbi = nbas(isym)
!  do ib=1,nbi
!    write(u6,*) 'orbital',isym,ib
!    write(u6,'(4ES19.12)') (CMO(imo+i),i=1,nbi)
!  imo = imo+nbi
!  end do
!end do
!----------------------------------------------------------------------*
!     Localize the inactive and virtual orbitals                       *
!----------------------------------------------------------------------*
call Cho_x_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,nAsh,nSsh,CMO,nCMO)
if (irc /= 0) then
  write(u6,*) 'Localization failed. The AFRE option cannot be used'
  call Abend()
end if
!write(u6,*) 'molecular orbitals after localization'
!imo = 0
!do isym=1,nsym
!  nbi = nbas(isym)
!  do ib=1,nbi
!    write(u6,*) 'orbital',isym,ib
!    write(u6,'(4ES19.12)') (CMO(imo+i),i=1,nbi)
!  imo = imo+nbi
!  end do
!end do
!----------------------------------------------------------------------*
!     Compute Mulliken atomic charges for each center and              *
!     each orbital.                                                    *
!----------------------------------------------------------------------*

nb2 = 0
do isym=1,nsym
  nb2 = nb2+nbas(isym)*(nbas(isym)+1)/2
end do
!write(u6,*) 'Starting the calculation',nb2
do i=1,nb2
  DPQ(i) = Zero
end do
ib = 0
imo0 = 0
ipq0 = 0
do isym=1,nsym
  nbi = nbas(isym)
  nfi = nfro(isym)
  nin = nish(isym)
  imo = imo0+nbi*nfi
  if (nin /= 0) then
    do i=1,nin
      labfro(i) = 0
    end do
    do ni=1,nin
      !write(u6,*) 'loop over sym and inactive orbitals',isym,ni
      ipq = ipq0
      ipq1 = 0
      do np=1,nbi
        do nq=1,np
          ipq = ipq+1
          ipq1 = ipq1+1
          DPQ(ipq1) = CMO(imo+np)*CMO(imo+nq)*SMAT(ipq)
        end do
      end do
      ! DPQ is the charge matrix for orbital ni in symmetry isym
      ! Now add non-diagonal elements to the diagonal
      ipq1 = 0
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        iqq = 0
        do nq=1,np
          iqq = iqq+nq
          ipq1 = ipq1+1
          if (np /= nq) then
            DPQ(ipp) = DPQ(ipp)+DPQ(ipq1)
            DPQ(iqq) = DPQ(iqq)+DPQ(ipq1)
          end if
        end do
      end do
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        !write(u6,*) 'diagonal element',ipp,DPQ(ipp)
      end do

      ! The diagonal now contains the charges for each basis function
      ! Add charges for basis functions centered on the selected atoms
      ! First check that the sum is equal to one
      chksum = Zero
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        chksum = chksum+DPQ(ipp)
      end do
      if (abs(chksum-One) > 1.0e-8_wp) then
        write(u6,*) 'Error on Checksum in Afreez. Value is not equal to 1:',isym,ni,chksum
        write(u6,*) 'Freezing extra orbitals in CASPT2 stops.'
        call Abend()
      end if
      ! Add diagonal elements that belong to selected atoms
      selch = Zero
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        do iname=1,lnfro
          if (name(ib+np)(1:4) == namfro(iname)) selch = selch+DPQ(ipp)
        end do
      end do
      if (abs(selch) < thrfr) labfro(ni) = 1
      !write(u6,*) selch
      imo = imo+nbi
    end do
    ! Sort the inactive CMO's such that frozen orbitals are first.
    nfro1 = nfro(isym)
    do ni=1,nin
      if (labfro(ni) == 1) then
        ! Exchange this orbital with the first inactive orbital
        ist1 = nfro(isym)*nbi+imo0
        ist2 = (nfro1+ni-1)*nbi+imo0
        !write(u6,*) 'nfro,nish',nfro(isym),nish(isym),ist1,ist2
        do np=1,nbi
          Swap = CMO(ist1+np)
          CMO(ist1+np) = CMO(ist2+np)
          CMO(ist2+np) = Swap
        end do
        nfro(isym) = nfro(isym)+1
        nish(isym) = nish(isym)-1
      end if
    end do
  end if
  ipq0 = ipq0+nbi*(nbi+1)/2
  imo0 = imo0+nbi**2
  ib = ib+nbi
end do
! Now sort virtual orbitals
! Orbitals with too low population on selected atoms will be deleted
do i=1,nb2
  DPQ(i) = Zero
end do
ib = 0
imo0 = 0
ipq0 = 0
do isym=1,nsym
  nbi = nbas(isym)
  ndi = ndel(isym)
  nsi = nssh(isym)
  nssh(isym) = 0
  ndel(isym) = ndi+nsi
  imo = imo0+nbi*(nfro(isym)+nish(isym)+nash(isym))
  if (nsi /= 0) then
    do i=1,nsi
      labfro(i) = 0
    end do
    do ni=1,nsi
      !write(u6,*) 'loop over sym and secondary orbitals',isym,ni
      ipq = ipq0
      ipq1 = 0
      do np=1,nbi
        do nq=1,np
          ipq = ipq+1
          ipq1 = ipq1+1
          DPQ(ipq1) = CMO(imo+np)*CMO(imo+nq)*SMAT(ipq)
        end do
      end do
      ! DPQ is the charge matrix for orbital ni in symmetry isym
      ! Now add non-diagonal elements to the diagonal
      ipq1 = 0
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        iqq = 0
        do nq=1,np
          iqq = iqq+nq
          ipq1 = ipq1+1
          if (np /= nq) then
            DPQ(ipp) = DPQ(ipp)+DPQ(ipq1)
            DPQ(iqq) = DPQ(iqq)+DPQ(ipq1)
          end if
        end do
      end do
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        !write(u6,*) 'diagonal element',ipp,DPQ(ipp)
      end do

      ! The diagonal now contains the charges for each basis function
      ! Add charges for basis functions centered on the selected atoms
      ! First check that the sum is equal to one
      chksum = Zero
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        chksum = chksum+DPQ(ipp)
      end do
      if (abs(chksum-One) > 1.0e-8_wp) then
        write(u6,*) 'Error on Checksum in Afreez. Value is not equal to 1:',isym,ni,chksum
        write(u6,*) 'Deleting extra orbitals in CASPT2 stops.'
        call Abend()
      end if
      !write(u6,*) 'Checksum',isym,ni,chksum
      ! Add diagonal elements that belong to selected atoms
      selch = Zero
      ipp = 0
      do np=1,nbi
        ipp = ipp+np
        do iname=1,lnfro
          if (name(ib+np)(1:4) == namfro(iname)) selch = selch+DPQ(ipp)
        end do
      end do
      if (abs(selch) > thrde) labfro(ni) = 1
      !write(u6,*) selch
      imo = imo+nbi
    end do
    ! Sort the CMO's such that secondary orbitals are first.
    do ni=1,nsi
      if (labfro(ni) == 1) then
        ! Exchange this orbital with the first deleted orbital
        ist1 = (nfro(isym)+nish(isym)+nash(isym)+nssh(isym))*nbi+imo0
        ist2 = (nfro(isym)+nish(isym)+nash(isym)+ni-1)*nbi+imo0
        do np=1,nbi
          Swap = CMO(ist1+np)
          CMO(ist1+np) = CMO(ist2+np)
          CMO(ist2+np) = Swap
        end do
        !write(u6,*) 'Orbital number',ni,ist1,ist2
        !write(u6,'(4ES19.12)') (CMO(ist1+np),np=1,nbi)
        !write(u6,'(4ES19.12)') (CMO(ist2+np),np=1,nbi)
        ndel(isym) = ndel(isym)-1
        nssh(isym) = nssh(isym)+1
      end if
    end do
  end if
  ipq0 = ipq0+nbi*(nbi+1)/2
  imo0 = imo0+nbi**2
  ib = ib+nbi
end do

!imo = 0
!do isym=1,nsym
!  nbi = nbas(isym)
!  do ib=1,nbi
!    write(u6,*) 'orbital',isym,ib
!    write(u6,'(4ES19.12)') (CMO(imo+i),i=1,nbi)
!  imo = imo+nbi
!  end do
!end do
! Write the resorted MO's back to JobIph

if (IFQCAN /= 0) IFQCAN = 0 ! MOs to be recanonicalized on exit

call MMA_DEALLOCATE(SMAT)

! Avoid unused argument warnings
if (.false.) call Unused_integer(NCMO)

end subroutine AFREEZ
