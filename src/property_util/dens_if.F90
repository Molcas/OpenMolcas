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

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: i_root
real(kind=wp), intent(_OUT_) :: CA(*), CB(*), OCCA(*), OCCB(*)
integer(kind=iwp) :: i, iA, iAC, iAC2, iad15, ii, IMO, IOCC, ip, ip1, ip2, ipAM1, ipAM2, ipC, ipDA, ipDB, ipDS, ipDT, ipUnity, &
                     ipVB, iS, J, nAct
real(kind=wp) :: Dum(1), OCCNO
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "casvb.fh"
#include "WrkSpc.fh"

call GetMem('DS','ALLO','REAL',ipDS,NACPAR)
call GetMem('DT','ALLO','REAL',ipDT,NACPAR)
call GetMem('DA','ALLO','REAL',ipDA,NACPAR)
call GetMem('DB','ALLO','REAL',ipDB,NACPAR)
call GetMem('CMO','ALLO','REAL',ipC,NTOT2)
call GetMem('UNITY','ALLO','REAL',ipUnity,NAC**2)

! READ IN ORBITALS
! Averaged...
if (iOrbTyp /= 2) iad15 = IADR15(2)
! Canonical...
if (iOrbTyp == 2) iad15 = IADR15(9)

! Read-in orbitals from Jobiph following instructions from previous lines...
call DDAFile(JOBIPH,2,Work(ipC),NTOT2,IAD15)

! COEN WANTED IT AS A BLOCKED MATRIX, SO HERE THEY COME...
ip1 = ipC
ip2 = 1
call dcopy_(nTOT**2,[Zero],0,CA,1)
call dcopy_(nTOT**2,[Zero],0,CB,1)
do iS=1,nSym
  do i=1,nbas(is)
    call dcopy_(nbas(is),Work(ip1),1,CA(ip2),1)
    call dcopy_(nbas(is),Work(ip1),1,CB(ip2),1)
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
    call dcopy_(nTOT,[Zero],0,OCCB,1)
  end if
  ! Canonical...
  if (iOrbTyp == 2) then
    ! SIMPLY ZERO ALPHA and BETA OCC NOS
    call dcopy_(nTOT,[Zero],0,OCCA,1)
    call dcopy_(nTOT,[Zero],0,OCCB,1)
  end if

  if (IFVB /= 0) then
    call GetMem('VB','ALLO','REAL',ipVB,NAC*NAC)
    call getvb2mo_cvb(Work(ipVB))
    call GetMem('ACTMO1','ALLO','REAL',ipAM1,NTOT*NAC)
    call GetMem('ACTMO2','ALLO','REAL',ipAM2,NTOT*NAC)
    ! Gather active MOs ...
    ! Also count no of active electrons ...
    iAC = 0
    IMO = 1
    IOCC = 1
    OCCNO = Zero
    nAct = 0
    do iS=1,nSym
      call dcopy_(nTOT*nash(is),CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1,Work(IAC+ipAM1),1)
      do J=0,NASH(IS)-1
        OCCNO = OCCNO+OCCA(J+NISH(IS)+NFRO(IS)+IOCC)
      end do
      nAct = nAct+NASH(iS)
      IAC = IAC+NASH(iS)*NTOT
      IMO = IMO+nBas(is)*ntot
      IOCC = IOCC+NBAS(iS)
    end do

    call DGEMM_('N','N',NTOT,NAC,NAC,One,Work(ipAM1),NTOT,Work(ipVB),NAC,Zero,Work(ipAM2),NTOT)

    ! Scatter active MOs ...
    ! Also reset active occ nos - we choose nel/nact since that
    ! gives as much meaning as anything ...

    iAC = 0
    IMO = 1
    IOCC = 1
    do iS=1,nSym
      call dcopy_(nTOT*nash(is),Work(IAC+ipAM2),1,CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1)
      call dcopy_(nash(is),[OCCNO/real(nAct,kind=wp)],0,OCCA(NISH(iS)+NFRO(IS)+IOCC),1)
      IAC = IAC+NASH(iS)*NTOT
      IMO = IMO+nBas(is)*ntot
      IOCC = IOCC+NBAS(iS)
    end do
    call GetMem('ACTMO1','FREE','REAL',ipAM1,NTOT*NAC)
    call GetMem('ACTMO2','FREE','REAL',ipAM2,NTOT*NAC)
    call GetMem('VB','FREE','REAL',ipVB,NAC*NAC)
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
    call DDAFile(JOBIPH,2,Work(ipDS),NACPAR,IAD15)
    call DDAFile(JOBIPH,2,Work(ipDT),NACPAR,IAD15)
    call DDAFile(JOBIPH,0,Dum,NACPR2,IAD15)
    call DDAFile(JOBIPH,0,Dum,NACPR2,IAD15)
  end do

  ! CREATE SPIN DENSITIES

  do i=0,NACPAR-1
    Work(ipDA+i) = Half*(Work(ipDS+i)+Work(ipDT+i))
    Work(ipDB+i) = Half*(Work(ipDS+i)-Work(ipDT+i))
  end do


  ! DIAGONALIZE THE SPIN DENSITIES

  ! FIRST ALPHA

  call dcopy_(nac**2,[Zero],0,Work(ipUnity),1)
  call dcopy_(nac,[One],0,Work(ipUnity),nac+1)
  call Jacob(Work(ipDA),Work(ipUNITY),NAC,NAC)

  ! TRANSFORM THE ACTIVE ORBITALS

  iAC = 0
  iAC2 = 1
  ip = ipC
  do iS=1,nSym
    if (nBas(iS) /= 0) then
      ip = ip+nBas(iS)*nIsh(iS)
      iAC2 = iAC2+nish(is)*nTot
      if (NASH(IS) > 0) call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),One,Work(ip),NBAS(IS),WORK(IPUNITY+IAC),NAC,Zero,CA(iAC2), &
                                    NTOT)
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
    call dcopy_(nBAS(is),[Zero],0,OCCA(ip),1)
    call dcopy_((nFro(is)+nish(is)),[One],0,OCCA(ip),1)
    ip = ip+nFro(is)+nish(is)
    do iA=1,nash(is)
      ii = ii+1
      i = i+ii
      OCCA(ip) = Work(ipDA+i-1)
      ip = ip+1
    end do
    ip = ip+nbas(is)-nFro(is)-nish(is)-nash(is)
  end do

  ! THEN ONCE AGAIN FOR BETA....

  call dcopy_(nac**2,[Zero],0,Work(ipUnity),1)
  call dcopy_(nac,[One],0,Work(ipUnity),nac+1)

  call Jacob(Work(ipDB),Work(ipUNITY),NAC,NAC)

  iAC = 0
  iAC2 = 1
  ip = ipC
  do iS=1,nSym
    if (nbas(is) /= 0) then
      ip = ip+nBas(is)*nIsh(is)
      iAC2 = iAC2+nish(is)*NTOT
      if (NASH(IS) > 0) call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),One,Work(ip),NBAS(IS),WORK(IPUNITY+IAC),NAC,Zero,CB(iac2), &
                                    NTOT)
      iAC = iAC+NASH(IS)*NAC+NASH(IS)
      iAC2 = iAC2+(nbas(is)-nish(is))*NTOT+nbas(is)
      ip = ip+nbas(is)*(nbas(is)-nish(is))
    end if
  end do

  ip = 1
  i = 0
  ii = 0
  do iS=1,nSYM
    call dcopy_(nBAS(is),[Zero],0,OCCB(ip),1)
    call dcopy_((nFro(is)+nish(is)),[One],0,OCCB(ip),1)
    ip = ip+nFro(is)+nish(is)
    do iA=1,nash(is)
      ii = ii+1
      i = i+ii
      OCCB(ip) = Work(ipDB+i-1)
      ip = ip+1
    end do
    ip = ip+nbas(is)-nFro(is)-nish(is)-nash(is)
  end do
end if

! OK, CLEAN UP

call GetMem('DS','FREE','REAL',ipDS,NACPAR)
call GetMem('DT','FREE','REAL',ipDT,NACPAR)
call GetMem('DA','FREE','REAL',ipDA,NACPAR)
call GetMem('DB','FREE','REAL',ipDB,NACPAR)
call GetMem('CMO','FREE','REAL',ipC,NTOT2)
call GetMem('UNITY','FREE','REAL',ipUnity,NAC**2)

return

end subroutine Dens_IF
