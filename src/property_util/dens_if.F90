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
! Copyright (C) 1999, Anders Bernhardsson                              *
!***********************************************************************

subroutine Dens_IF(i_root,CA,CB,OCCA,OCCB)
! A small stupid interface for creating the alpha and beta
! occupation numbers and corresponding molecular orbitals,
! used in MOLDEN.
!
! EAW 990118

use casvb_global, only: ifvb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: i_root
real(kind=wp), intent(_OUT_) :: CA(*), CB(*), OCCA(*), OCCB(*)
integer(kind=iwp) :: i, iA, iAC, iAC2, iad15, ii, IMO, IOCC, ip, ip1, ip2, iS, J, nAct
real(kind=wp) :: Dum(1), OCCNO
real(kind=wp), allocatable :: AM1(:,:), AM2(:,:), C(:), DA(:), DB(:), DS(:), DT(:), Unity(:), VB(:,:)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

call mma_allocate(DS,NACPAR,label='DS')
call mma_allocate(DT,NACPAR,label='DT')
call mma_allocate(DA,NACPAR,label='DA')
call mma_allocate(DB,NACPAR,label='DB')
call mma_allocate(C,NTOT2,label='C')
call mma_allocate(Unity,NAC*NAC,label='UNITY')

! READ IN ORBITALS
! Averaged...
if (iOrbTyp /= 2) iad15 = IADR15(2)
! Canonical...
if (iOrbTyp == 2) iad15 = IADR15(9)

! Read-in orbitals from Jobiph following instructions from previous lines...
call DDAFile(JOBIPH,2,C,NTOT2,IAD15)

! COEN WANTED IT AS A BLOCKED MATRIX, SO HERE THEY COME...
ip1 = 1
ip2 = 1
CA(:nTOT**2) = Zero
CB(:nTOT**2) = Zero
do iS=1,nSym
  do i=1,nbas(is)
    call dcopy_(nbas(is),C(ip1),1,CA(ip2),1)
    call dcopy_(nbas(is),C(ip1),1,CB(ip2),1)
    ip1 = ip1+nbas(is)
    ip2 = ip2+NTOT
  end do
  ip2 = ip2+nbas(is)
end do

! --- FOLLOWING FOR STANDARD (PSEUDO/AVERAGED/CANONICAL) ORBITALS ---
if (i_root == 0) then
  ! Averaged...
  if (iOrbTyp /= 2) then
    ! SIMPLY READ OCC NOS AS ALPHAS AND ZERO BETA OCC NOS
    iad15 = IADR15(2)
    Dum(1) = Zero
    call DDAFile(JOBIPH,0,Dum,NTOT2,IAD15)
    call DDAFILE(JOBIPH,2,OCCA,NTOT,IAD15)
    OCCB(:nTOT) = Zero
  end if
  ! Canonical...
  if (iOrbTyp == 2) then
    ! SIMPLY ZERO ALPHA and BETA OCC NOS
    OCCA(:nTOT) = Zero
    OCCB(:nTOT) = Zero
  end if

  if (IFVB /= 0) then
    call mma_allocate(VB,NAC,NAC,label='VB')
    call getvb2mo_cvb(VB)
    call mma_allocate(AM1,NTOT,NAC,label='ACTMO1')
    call mma_allocate(AM2,NTOT,NAC,label='ACTMO2')
    ! Gather active MOs ...
    ! Also count no of active electrons ...
    iAC = 1
    IMO = 1
    IOCC = 1
    OCCNO = Zero
    nAct = 0
    do iS=1,nSym
      call dcopy_(nTOT*nash(is),CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1,AM1(:,IAC:IAC+NASH(iS)-1),1)
      do J=0,NASH(IS)-1
        OCCNO = OCCNO+OCCA(J+NISH(IS)+NFRO(IS)+IOCC)
      end do
      nAct = nAct+NASH(iS)
      IAC = IAC+NASH(iS)
      IMO = IMO+nBas(is)*ntot
      IOCC = IOCC+NBAS(iS)
    end do

    call DGEMM_('N','N',NTOT,NAC,NAC,One,AM1,NTOT,VB,NAC,Zero,AM2,NTOT)

    ! Scatter active MOs ...
    ! Also reset active occ nos - we choose nel/nact since that
    ! gives as much meaning as anything ...

    iAC = 1
    IMO = 1
    IOCC = 1
    do iS=1,nSym
      call dcopy_(nTOT*nash(is),AM2(:,IAC:IAC+NASH(iS)-1),1,CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1)
      call dcopy_(nash(is),[OCCNO/real(nAct,kind=wp)],0,OCCA(NISH(iS)+NFRO(IS)+IOCC),1)
      IAC = IAC+NASH(iS)
      IMO = IMO+nBas(is)*ntot
      IOCC = IOCC+NBAS(iS)
    end do
    call mma_deallocate(VB)
    call mma_deallocate(AM1)
    call mma_deallocate(AM2)
  end if

else
!                                                                      *
!***********************************************************************
!                                                                      *
! --- FOLLOWING FOR NATURAL *SPIN* ORBITALS ---
!
! READ IN DENSITIES

  iad15 = IADR15(3)
  Dum(1) = Zero
  do i=1,i_root
    call DDAFile(JOBIPH,2,DS,NACPAR,IAD15)
    call DDAFile(JOBIPH,2,DT,NACPAR,IAD15)
    call DDAFile(JOBIPH,0,Dum,NACPR2,IAD15)
    call DDAFile(JOBIPH,0,Dum,NACPR2,IAD15)
  end do

  ! CREATE SPIN DENSITIES

  DA(:) = Half*(DS(:)+DT(:))
  DB(:) = Half*(DS(:)-DT(:))

  ! DIAGONALIZE THE SPIN DENSITIES

  ! FIRST ALPHA

  call unitmat(Unity,nac)
  call Jacob(DA,UNITY,NAC,NAC)

  ! TRANSFORM THE ACTIVE ORBITALS

  iAC = 1
  iAC2 = 1
  ip = 1
  do iS=1,nSym
    if (nBas(iS) /= 0) then
      ip = ip+nBas(iS)*nIsh(iS)
      iAC2 = iAC2+nish(is)*nTot
      if (NASH(IS) > 0) call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),One,C(ip),NBAS(IS),UNITY(IAC),NAC,Zero,CA(iAC2),NTOT)
      iAC = iAC+NASH(IS)*NAC+NASH(IS)
      iAC2 = iAC2+(nbas(is)-nish(is))*NTOT+nbas(is)
      ip = ip+nbas(is)*(nbas(is)-nish(is))
    end if
  end do

  ! COPY OCCUPATION NUMBERS

  ip = 1
  i = 0
  ii = 0
  do iS=1,nSYM
    OCCA(ip:ip+nBAS(is)-1) = Zero
    OCCA(ip:ip+nFro(is)+nish(is)-1) = One
    ip = ip+nFro(is)+nish(is)
    do iA=1,nash(is)
      ii = ii+1
      i = i+ii
      OCCA(ip) = DA(i)
      ip = ip+1
    end do
    ip = ip+nbas(is)-nFro(is)-nish(is)-nash(is)
  end do

  ! THEN ONCE AGAIN FOR BETA....

  call unitmat(Unity,nac)

  call Jacob(DB,UNITY,NAC,NAC)

  iAC = 1
  iAC2 = 1
  ip = 1
  do iS=1,nSym
    if (nbas(is) /= 0) then
      ip = ip+nBas(is)*nIsh(is)
      iAC2 = iAC2+nish(is)*NTOT
      if (NASH(IS) > 0) call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),One,C(ip),NBAS(IS),UNITY(IAC),NAC,Zero,CB(iac2),NTOT)
      iAC = iAC+NASH(IS)*NAC+NASH(IS)
      iAC2 = iAC2+(nbas(is)-nish(is))*NTOT+nbas(is)
      ip = ip+nbas(is)*(nbas(is)-nish(is))
    end if
  end do

  ip = 1
  i = 0
  ii = 0
  do iS=1,nSYM
    OCCB(ip:ip+nBAS(is)-1) = Zero
    OCCB(ip:ip+nFro(is)+nish(is)-1) = One
    ip = ip+nFro(is)+nish(is)
    do iA=1,nash(is)
      ii = ii+1
      i = i+ii
      OCCB(ip) = DB(i)
      ip = ip+1
    end do
    ip = ip+nbas(is)-nFro(is)-nish(is)-nash(is)
  end do
end if

! OK, CLEAN UP

call mma_deallocate(DS)
call mma_deallocate(DT)
call mma_deallocate(DA)
call mma_deallocate(DB)
call mma_deallocate(C)
call mma_deallocate(Unity)

return

end subroutine Dens_IF
