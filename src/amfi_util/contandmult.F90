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

subroutine contandmult(Lhigh,AIMP,oneonly,numballcart,LUPROP,ifinite,onecart,onecontr,oneoverR3,iCenter)

use AMFI_global, only: charge, contrarray, icore, ikeeplist, ikeeporb, incrLM, ipowxyz, iredLM, iredoffunctnew, itotalperIR, &
                       Lmax, Loffunction, Moffunction, MxcontL, MxprimL, ncontrac, ncontrac_keep, nprimit, nrtofiperIR, numbofsym, &
                       shiftIRED, shiftIRIR
use index_functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lhigh, numballcart, LUPROP, ifinite, iCenter
logical(kind=iwp), intent(in) :: AIMP, oneonly
real(kind=wp), intent(inout) :: onecart(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,3)
real(kind=wp), intent(out) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3,Lmax)
real(kind=wp), intent(in) :: oneoverR3(MxprimL*(MxprimL+1)/2,Lmax)
integer(kind=iwp) :: I, icartfirst, icartsec, ind1, ind2, ipntnew, ipntold, ipowx, ipowy, ipowz, ired1, ired2, iredfirst, &
                     iredired, iredsec, irun, irun1, irun2, L, length3, Lrun, Mfirst, mrun, Msec, norb1, norb2, norbsh1, norbsh2
character(len=8) :: xa(4), ya(4), za(4)
real(kind=wp), allocatable :: Dummy(:), OCA(:,:), OCA2(:,:), OCA3(:,:)

!bs get back the real number of functions for the finite nucleus
if (ifinite == 2) ncontrac(0) = ncontrac_keep
!#######################################################################
!bs subroutine to contract radial one-electron integrals
!bs and multiply them with angular factors
!#######################################################################
xa(:) = ['********','        ','ANTISYMM','X1SPNORB']
ya(:) = ['********','        ','ANTISYMM','Y1SPNORB']
za(:) = ['********','        ','ANTISYMM','Z1SPNORB']

!bs clean the arrays for cartesian integrals

length3 = numbalLcart*(numbalLcart+1)/2
call mma_allocate(OCA,Length3,3,Label='OCA')
call mma_allocate(OCA2,Length3,3,Label='OCA2')
call mma_allocate(Dummy,MxContL**2,Label='Dummy')
Dummy(:) = Zero
OCA(:,:) = Zero
OCA2(:,:) = Zero

!bs one-electron-integrals:
!bs 1. index: number of first contracted function
!bs 2. index: number of second contracted function
!bs 3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
!bs 4. index: L-value
!bs onecart(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,1),
!bs onecart(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,2),
!bs onecart(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,3)

!bs generate one-electron integrals for all L greater/equal 1
if (ifinite == 2) charge = Zero ! nuclear integrals are modelled for finite nucleus somewhere else
do L=1,Lhigh
  call contone(L,oneoverr3(:,L),onecontr(:,:,-Lmax,1,L),Lmax,contrarray(:,3,L),nprimit(L),ncontrac(L),MxcontL,Dummy, &
               onecart(:,:,:,L,1),onecart(:,:,:,L,2),onecart(:,:,:,L,3),charge,oneonly)
end do

!bs ********************************************************************
!bs now move all integrals to one big arrays for X,Y,Z
!bs ********************************************************************
do Lrun=1,Lhigh  !loop over L-values (integrals are diagonal in L)
  mrun = 0
  do Msec=-Lrun,Lrun    ! cartesian M-values  (Mfirst,Msec) with
    do Mfirst=-Lrun,Msec  ! Mfirst <= Msec (actually '=' does never appear as there is no L-component  in Ag

      !bs determine  if L_X L_Y or L_Z
      ipowx = ipowxyz(1,mfirst,Lrun)+ipowxyz(1,msec,Lrun)
      ipowy = ipowxyz(2,mfirst,Lrun)+ipowxyz(2,msec,Lrun)
      ipowz = ipowxyz(3,mfirst,Lrun)+ipowxyz(3,msec,Lrun)

      mrun = mrun+1
      !bs now determine the irreducible representations
      iredfirst = iredLM(Mfirst,Lrun)
      iredsec = iredLM(Msec,Lrun)
      !bs check out which IR is the lower one.
      if (iredfirst <= iredsec) then

        !bs calculate shift to get to the beginning of the block
        iredired = shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)+incrlm(Mfirst,Lrun)*itotalperIR(iredsec)+incrLM(Msec,Lrun)
        if ((mod(ipowx,2) == 0) .and. (mod(ipowy,2) == 1) .and. (mod(ipowz,2) == 1)) then
          do icartfirst=1,ncontrac(Lrun) ! loop first index
            do icartsec=1,ncontrac(Lrun) ! loop second index
              oca(iredired+icartsec,1) = oca(iredired+icartsec,1)+onecart(icartfirst,icartsec,mrun,Lrun,1)
            end do
            !bs  shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredsec)
          end do
        end if
        if ((mod(ipowx,2) == 1) .and. (mod(ipowy,2) == 0) .and. (mod(ipowz,2) == 1)) then
          do icartfirst=1,ncontrac(Lrun) ! loop first index
            do icartsec=1,ncontrac(Lrun) ! loop second index
              oca(iredired+icartsec,2) = oca(iredired+icartsec,2)+onecart(icartfirst,icartsec,mrun,Lrun,2)
            end do
            !bs shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredsec)
          end do
        end if
        if ((mod(ipowx,2) == 1) .and. (mod(ipowy,2) == 1) .and. (mod(ipowz,2) == 0)) then
          do icartfirst=1,ncontrac(Lrun) ! loop first index
            do icartsec=1,ncontrac(Lrun) ! loop second index
              oca(iredired+icartsec,3) = oca(iredired+icartsec,3)+onecart(icartfirst,icartsec,mrun,Lrun,3)
            end do
            !bs shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredsec)
          end do
        end if
      else if (iredfirst > iredsec) then
        !bs In this case, indices are exchanged with respect to former
        !bs symmetry of blocks. Therefore, there will be a minus sign

        !bs calculate shift to get to the beginning of the block
        iredired = shiftIRIR((iredfirst*iredfirst-iredfirst)/2+iredsec)+incrLM(Msec,Lrun)*itotalperIR(iredfirst)+incrLM(Mfirst,Lrun)
        if ((mod(ipowx,2) == 0) .and. (mod(ipowy,2) == 1) .and. (mod(ipowz,2) == 1)) then
          do icartsec=1,ncontrac(Lrun)     !loop second index
            do icartfirst=1,ncontrac(Lrun) !loop first index
              oca(iredired+icartfirst,1) = oca(iredired+icartfirst,1)-onecart(icartsec,icartfirst,mrun,Lrun,1)
            end do
            !bs shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredfirst)
          end do
        end if
        if ((mod(ipowx,2) == 1) .and. (mod(ipowy,2) == 0) .and. (mod(ipowz,2) == 1)) then
          do icartsec=1,ncontrac(Lrun)     !loop second index
            do icartfirst=1,ncontrac(Lrun) !loop first index
              oca(iredired+icartfirst,2) = oca(iredired+icartfirst,2)-onecart(icartsec,icartfirst,mrun,Lrun,2)
            end do
            !bs shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredfirst)
          end do
        end if
        if ((mod(ipowx,2) == 1) .and. (mod(ipowy,2) == 1) .and. (mod(ipowz,2) == 0)) then
          do icartsec=1,ncontrac(Lrun)     !loop  second index
            do icartfirst=1,ncontrac(Lrun) !loop first index
              oca(iredired+icartfirst,3) = oca(iredired+icartfirst,3)-onecart(icartsec,icartfirst,mrun,Lrun,3)
            end do
            !bs shift pointer by number of functions in IR
            iredired = iredired+itotalperIR(iredfirst)
          end do
        end if
      end if
    end do
  end do
end do

!bs copy integrals on arrays with no symmetry blocking at all
!bs which means huge triangular matrices
irun = 0
do norb2=1,numballcart
  ired2 = iredoffunctnew(norb2)
  norbsh2 = norb2-shiftIRED(ired2)
  do norb1=1,norb2
    ired1 = iredoffunctnew(norb1)
    norbsh1 = noRb1-shiftIRED(ired1)
    irun = irun+1
    iredired = shiftIRIR((ired2*ired2-ired2)/2+ired1)
    if (ired1 /= ired2) then
      oca2(irun,:) = oca(iredired+norbsh2+(norbsH1-1)*itotalperIR(IREd2),:)
    else
      oca2(irun,:) = oca(iredired+norbsh2*(norbsH2-1)/2+norbsh1,:)
    end if
  end do
end do
if (.not. AIMP) then
  ! write a hermit-like file   b.s. 4.10.96
  !BS write(u6,*) 'number of orbitals ',numbalLcart
  !BS write(u6,*) 'length of triangular matrix ', length3
  !BS This was removed and will be done in SEWARD
  !BS open(LUPROP,status='UNKNOWN',form='UNFORMATTED',file='AOPROPER_MF')
  !BS rewind(LUPROP)
  write(LUPROP) iCenter
  write(LUPROP) xa,numbofsym,(nrtofiperIR(I),i=1,numbofsym),numballcart,(Loffunction(I),I=1,numballcart), &
                (Moffunction(I),I=1,numballcart),Lhigh,(ncontrac(I),I=0,Lhigh)
  write(LUPROP) (oca2(irun,1),irun=1,length3)
  write(LUPROP) Ya
  write(LUPROP) (oca2(irun,2),irun=1,length3)
  write(LUPROP) Za
  write(LUPROP) (oca2(irun,3),irun=1,length3)
  !BS close(LUPROP)
else
  !bs reorder for AIMP
  !bs write(u6,*) 'reorder integrals for AIMP'
  length3 = ikeeporb*(ikeeporb+1)/2
  call mma_allocate(OCA3,length3,3,Label='OCA3')
  OCA3(:,:) = Zero
  !bs write(u6,*) 'number of orbitals ',ikeeporb
  !bs write(u6,*) 'length of triangular matrix ', length3
  do irun2=1,ikeeporb
    do irun1=1,irun2
      ind2 = ikeeplist(irun2)
      ind1 = ikeeplist(irun1)
      ipntold = iTri(ind1,ind2)
      ipntnew = iTri(irun1,irun2)

      oca3(ipntnew,:) = oca2(ipntold,:)
    end do
  end do
  !BS write(u6,*) 'transfered to new blocks'
  !BS LUPROP = 19
  !BS open(LUPROP,status='UNKNOWN',form='UNFORMATTED',file='AOPROPER_MF')
  !BS rewind(LUPROP)

  write(LUPROP) iCenter
  write(LUPROP) xa,numbofsym,(nrtofiperIR(I),i=1,numbofsym),ikeeporb,(Loffunction(ikeeplist(i)),i=1,ikeeporb), &
                (Moffunction(ikeeplist(i)),I=1,ikeeporb),Lhigh,((ncontrac(I)-icore(I)),I=0,Lhigh)
  write(LUPROP) (oca3(irun,1),irun=1,length3)
  write(LUPROP) ya
  write(LUPROP) (oca3(irun,2),irun=1,length3)
  write(LUPROP) za
  write(LUPROP) (oca3(irun,3),irun=1,length3)

  call mma_deallocate(OCA3)
  !BS close(LUPROP)
end if

!bs that is it!!

call mma_deallocate(OCA2)
call mma_deallocate(OCA)
call mma_deallocate(Dummy)

return

end subroutine contandmult
