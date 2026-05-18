!ERI0***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine wfnsizes()
!***********************************************************************
!
! Compute various orbital sizes
!
!***********************************************************************

use Molcas, only: MxAct, MxIna, MxOrb
use caspt2_global, only: NTAT, NTORB, NPREF, NDREF
use caspt2_module, only: iSCF, iSpin, MxExt, nActEl, nAmx, nAshT, nBasT, nBMx, nBSqT, nBTri, nIMx, nInaBx, nIshT, nOMx, nOrbT, &
                         nOSqT, nOTri, nRas1T, nRas2T, nRas3T, nSecBx, nSMx, nSshT, nSym, nIes, nAes, nSes, nOsh, nFroT, nAsh, &
                         nSsh, nDel, nOrb, nIsh, nFro, nRas1, nRas2, nRas3, nBas, OrbNam, IINAIS, iExtIS, iiSym, iaSym, ISNAM
use caspt2_module, only: nG1, nG2, nG3Tot
use definitions, only: iwp, u6

implicit none
integer(kind=iwp) NASHT2
integer(kind=iwp) NI, NR1, NR2, NR3, NS, N123
integer(kind=iwp) I, iSym
! Table sizes
integer(kind=iwp) IIABS, ITABS
integer(kind=iwp) IS, IO
integer(kind=iwp) ITOT, IINA, IEXT

NFROT = 0
NISHT = 0
NASHT = 0
NRAS1T = 0
NRAS2T = 0
NRAS3T = 0
NSSHT = 0
NORBT = 0
NBAST = 0
NOSQT = 0
NBSQT = 0
NIMX = 0
NAMX = 0
NSMX = 0
NOMX = 0
NBMX = 0
do ISYM=1,NSYM
  NIES(ISYM) = NISHT
  NAES(ISYM) = NASHT
  NSES(ISYM) = NSSHT
  NOSH(ISYM) = NISH(ISYM)+NASH(ISYM)
  NSSH(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NOSH(ISYM)-NDEL(ISYM)
  NORB(ISYM) = NOSH(ISYM)+NSSH(ISYM)
  NORBT = NORBT+NORB(ISYM)
  NOSQT = NOSQT+NORB(ISYM)**2
  NBSQT = NBSQT+NBAS(ISYM)**2
  NFROT = NFROT+NFRO(ISYM)
  NISHT = NISHT+NISH(ISYM)
  NASHT = NASHT+NASH(ISYM)
  NRAS1T = NRAS1T+NRAS1(ISYM)
  NRAS2T = NRAS2T+NRAS2(ISYM)
  NRAS3T = NRAS3T+NRAS3(ISYM)
  NSSHT = NSSHT+NSSH(ISYM)
  NBAST = NBAST+NBAS(ISYM)
  NIMX = max(NIMX,NISH(ISYM))
  NAMX = max(NAMX,NASH(ISYM))
  NSMX = max(NSMX,NSSH(ISYM))
  NOMX = max(NOMX,NORB(ISYM))
  NBMX = max(NBMX,NBAS(ISYM))
end do
! Set RHS Boxes to maximum size
NINABX = NIMX
NSECBX = NSMX
NBTRI = (NBSQT+NBAST)/2
NOTRI = (NOSQT+NORBT)/2
! Size of orbital transformation arrays:
NTORB = 0
NTAT = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NR1 = NRAS1(ISYM)
  NR2 = NRAS2(ISYM)
  NR3 = NRAS3(ISYM)
  NS = NSSH(ISYM)
  N123 = NR1**2+NR2**2+NR3**2
  NTAT = NTAT+N123
  NTORB = NTORB+NI**2+N123+NS**2
end do
! Sizes of DREF and PREF arrays:
NDREF = 1
NPREF = 1
! Sizes of GAMMA1, GAMMA2, and GAMMA3:
NG1 = 1
NG2 = 1
NG3TOT = 1
if (NASHT > 0) then
  NDREF = (NASHT**2+NASHT)/2
  NASHT2 = NASHT**2
  NPREF = (NASHT2**2+NASHT2)/2
  NG1 = NASHT**2
  NG2 = NG1**2
  NG3TOT = ((NG1+2)*(NG1+1)*NG1)/6
end if

! Identify the wave function type
ISCF = 0
if ((ISPIN == NACTEL+1) .and. (NACTEL == NASHT)) ISCF = 2
if (NASHT == 0) ISCF = 1
if (NACTEL == 2*NASHT) ISCF = 1

!***********************************************************************
!
! Create orbital name vector
!
!***********************************************************************
IS = 0
ITOT = 0
IINA = 0
IEXT = 0
do ISYM=1,NSYM
  IO = 0
  do I=1,NFRO(ISYM)
    ITOT = ITOT+1
    IO = IO+1
    write(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)') 'Fr',ISYM,'.',IO
  end do
  do I=1,NISH(ISYM)
    ITOT = ITOT+1
    IINA = IINA+1
    IINAIS(IINA) = ITOT
    IO = IO+1
    write(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)') 'In',ISYM,'.',IO
  end do
  do I=1,NASH(ISYM)
    ITOT = ITOT+1
    IO = IO+1
    write(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)') 'Ac',ISYM,'.',IO
  end do
  do I=1,NSSH(ISYM)
    ITOT = ITOT+1
    IEXT = IEXT+1
    IEXTIS(IEXT) = ITOT
    IO = IO+1
    write(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)') 'Se',ISYM,'.',IO
    IS = IS+1
    ISNAM(IS) = ORBNAM(ITOT)
  end do
  do I=1,NDEL(ISYM)
    ITOT = ITOT+1
    IO = IO+1
    write(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)') 'De',ISYM,'.',IO
  end do
end do

!***********************************************************************
!
! Precompute table sizes
!
!***********************************************************************
IIABS = 0
ITABS = 0
do ISYM=1,NSYM
  do I=1,NORB(ISYM)
    if (I <= NISH(ISYM)) then
      IIABS = IIABS+1
      IISYM(IIABS) = ISYM
    else if (I <= NISH(ISYM)+NASH(ISYM)) then
      ITABS = ITABS+1
      IASYM(ITABS) = ISYM
    end if
  end do
end do

! Check consistency of the orbitals
if (NISHT > MXINA) then
  call WarningMessage(2,'Too many inactive orbitals.')
  write(u6,'(a,2i8)') ' NISHT >  MXINA:',NISHT,MXINA
  call Quit_OnUserError()
end if
if (NASHT > MXACT) then
  call WarningMessage(2,'Too many active orbitals.')
  write(u6,'(a,2i8)') ' NASHT > MXACT:',NASHT,MXACT
  call Quit_OnUserError()
end if
if (NSSHT > MXEXT) then
  call WarningMessage(2,'Too many secondary orbitals.')
  write(u6,'(a,2i8)') ' NSSHT > MXEXT:',NSSHT,MXEXT
  call Quit_OnUserError()
end if
if (NBAST > MXORB) then
  call WarningMessage(2,'Too many basis functions.')
  write(u6,'(a,2i8)') ' NBAST > MXORB:',NBAST,MXORB
  call Quit_OnUserError()
end if

! GG-Nov04  The following informations must be passed to the Cholesky
! transformation section through RunFile. COMMON blocks cannot be
! used due to several conflicts.
call Put_iScalar('nSym',nSym)
call Put_iArray('nFroPT',nFro,nSym)
call Put_iArray('nIsh',nIsh,nSym)
call Put_iArray('nAsh',nAsh,nSym)
call Put_iArray('nDelPT',nDel,nSym)
call Put_iArray('nBas',nBas,nSym)

end subroutine wfnsizes
