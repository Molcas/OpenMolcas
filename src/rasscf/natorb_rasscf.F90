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

subroutine NatOrb_RASSCF(CMOO,SCR1,SCR2,SMAT,CMON,OCCN)
! RASSCF program: version IBM-3090: Output section
!
! PURPOSE: Calculation of natural orbitals from the
!          density matrix. These orbitals are written onto JOBIPH
!          in position IADR15(12), followed by the occupation
!          numbers in IADR15(13).
!          The calculation is performed for each of the NROOT
!          density matrices obtained in an average CASSCF calc.
!          Called from MAIN before OUTCTL
!
!      ****** IBM 3090 MOLCAS Release: 90 02 22 ******

use Index_Functions, only: nTri_Elem
use rasscf_global, only: iADR15, KSDFT, lRoots, NACPAR, NACPR2
use SplitCas_Data, only: DoSPlitCas, lRootSplit
use PrintLevel, only: DEBUG, USUAL
use output_ras, only: IPRLOC
use general_data, only: JOBIPH, NASH, NBAS, NFRO, NISH, NSYM, NTOT, NTOT2
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMOO(*), SCR1(*), SCR2(*), SMAT(*), CMON(*), OCCN(*)
integer(kind=iwp) :: I, IA, IB, IBAS, ID, iDisk, IEND, II, IO, iOff, iPrLev, IST, iSTMO, ISTMO1, ISYM, J, JA, jDisk, jOff, kRoot, &
                     NA1, NAO, NB, NB2, NFI

IPRLEV = IPRLOC(7)
iDisk = IADR15(12)
jDisk = IADR15(3)

if (.not. DoSplitCAS) then
  do kRoot=1,lRoots
    if ((KSDFT == 'SCF') .and. (IPRLEV >= USUAL)) then
      write(u6,*)
      write(u6,'(6X,A,I3)') 'Natural orbitals and occupation numbers for root',kRoot
    end if
    call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
    call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
    call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
    call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
    call DBLOCK(SCR1)

    OCCN(1:NTOT) = Zero
    CMON(1:NTOT2) = CMOO(1:NTOT2)

    ID = 0
    ISTMO1 = 0
    IB = 0
    do ISYM=1,NSYM
      NB = NBAS(ISYM)
      NFI = NFRO(ISYM)+NISH(ISYM)
      NAO = NASH(ISYM)
      NB2 = NB**2
      NA1 = nTri_Elem(NAO)
      IO = IB+NFI
      ISTMO = ISTMO1+NB*NFI

      ! set occupation number of frozen and inactive orbitals

      OCCN(IB+1:IB+NFI) = Two

      ! Diagonalize the density matrix and transform orbitals

      if (NAO > 0) then
        call unitmat(SCR2,NAO)
        call JACOB(SCR1(ID+1),SCR2,NAO,NAO)
        II = 0
        do I=1,NAO
          II = II+I
          OCCN(IO+I) = SCR1(II+ID)
        end do
        IST = IO+1
        IEND = IO+NAO
        if ((KSDFT == 'SCF') .and. (IPRLEV >= USUAL)) write(u6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))') 'sym',iSym,':', &
                                                                                                        (OCCN(I),I=IST,IEND)
        call DGEMM_('N','N',NB,NAO,NAO,One,CMOO(ISTMO+1),NB,SCR2,NAO,Zero,CMON(ISTMO+1),NB)
      end if

      ID = ID+NA1
      ISTMO1 = ISTMO1+NB2
      IB = IB+NB
    end do

    ! ORTHOGONALIZE NEW MO'S

    call SUPSCH(SMAT,CMOO,CMON)
    call ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)

    ! Place NOs wrt decreasing order of OCCN
    !***  GLM commented it off for this reordering is not consistent
    !GLM throughout Molcas
    !
    !iOff = 0
    !jOff = 0
    !do iSym=1,nSym
    !  NB = nBas(iSym)
    !  NAO = nAsh(iSym)
    !  IA = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    !  JA = 1+jOff+nFro(iSym)+nIsh(iSym)
    !  call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
    !  iOff = iOff+NB**2
    !  jOff = jOff+NB
    !end do

    if (IPRLEV >= DEBUG) then
      write(u6,*)
      write(u6,*) ' CMON in NATORB_RASSCF after ORDER_ARRAYS'
      write(u6,*) ' ---------------------'
      write(u6,*)
      ioff = 0
      do iSym=1,nSym
        iBas = nBas(iSym)
        if (iBas /= 0) then
          write(u6,*) 'Sym =',iSym
          do i=1,iBas
            write(u6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
          end do
          iOff = iOff+(iBas*iBas)
        end if
      end do
    end if

    ! Write new molecular orbitals and occupation numbers to JOBIPH

    call DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
    call DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)
  end do

else ! if DoSplitCAS (GLMJ)...
  if ((KSDFT == 'SCF') .and. (IPRLEV >= USUAL)) then
    write(u6,*)
    write(u6,'(6X,A,I3)') 'Natural orbitals and occupation numbers for root',lRootSplit
  end if
  call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
  call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
  call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
  call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
  call DBLOCK(SCR1)

  OCCN(1:NTOT) = Zero
  CMON(1:NTOT2) = CMOO(1:NTOT2)

  ID = 0
  ISTMO1 = 0
  IB = 0
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    NFI = NFRO(ISYM)+NISH(ISYM)
    NAO = NASH(ISYM)
    NB2 = NB**2
    NA1 = nTri_Elem(NAO)
    IO = IB+NFI
    ISTMO = ISTMO1+NB*NFI

    ! set occupation number of frozen and inactive orbitals

    OCCN(IB+1:IB+NFI) = Two

    ! Diagonalize the density matrix and transform orbitals

    if (NAO > 0) then
      call unitmat(SCR2,NAO)
      call JACOB(SCR1(ID+1),SCR2,NAO,NAO)
      II = 0
      do I=1,NAO
        II = II+I
        OCCN(IO+I) = SCR1(II+ID)
      end do
      IST = IO+1
      IEND = IO+NAO
      if ((KSDFT == 'SCF') .and. (IPRLEV >= USUAL)) write(u6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))') 'sym',iSym,':', &
                                                                                                      (OCCN(I),I=IST,IEND)
      call DGEMM_('N','N',NB,NAO,NAO,One,CMOO(ISTMO+1),NB,SCR2,NAO,Zero,CMON(ISTMO+1),NB)
    end if

    ID = ID+NA1
    ISTMO1 = ISTMO1+NB2
    IB = IB+NB
  end do

  ! ORTHOGONALIZE NEW MO'S

  call SUPSCH(SMAT,CMOO,CMON)
  call ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)

  ! Place NOs wrt decreasing order of OCCN

  iOff = 0
  jOff = 0
  do iSym=1,nSym
    NB = nBas(iSym)
    NAO = nAsh(iSym)
    IA = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    JA = 1+jOff+nFro(iSym)+nIsh(iSym)
    call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
    iOff = iOff+NB**2
    jOff = jOff+NB
  end do

  ! Write new molecular orbitals and occupation numbers to JOBIPH

  call DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
  call DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)

end if

end subroutine NatOrb_RASSCF
