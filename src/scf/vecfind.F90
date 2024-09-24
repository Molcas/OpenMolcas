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
! Copyright (C) Per-Olof Widmark                                       *
!               2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine VecFind(OccSet,FermSet,CharSet,SpinSet)
!***********************************************************************
!                                                                      *
! This routine figure out which set of starting orbitals are used.     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!          R. Lindh                                                    *
!          Uppsala University, Sweden, Feb 2017                        *
!          Remove Work                                                 *
!                                                                      *
!***********************************************************************
!                                                                      *
! Input vector mappings                                                *
!                                                                      *
! LstVec              InVec                                            *
! -1 Die               0 Core                                          *
!  0 Old SCF           1 NDDO (not used)                               *
!  1 Guessorb          2 Lumorb                                        *
!  2 Lumorb            3 Old density                                   *
!  3 Old Density       4 Restart (not used)                            *
!  4 Core              8 Old SCF                                       *
!  5 NDDO (use not)    9 Guessorb                                      *
!                                                                      *
!***********************************************************************

#ifdef _HDF5_
use mh5, only: mh5_fetch_attr
use InfSCF, only: FileOrb_ID
#endif
use InfSCF, only: iAu_ab, InVec, isHDF5, nD, nSym, nStOpt, SCF_FileOrb, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, nBas, LstVec, &
                  nOcc, nAufb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
logical OccSet
logical FermSet
logical CharSet
logical SpinSet
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
character*80 cList
logical Found
integer mynSym
integer mynBas(8)
integer mynOrb(8)
character*10 infoLbl
real*8, dimension(:), allocatable :: EOrb
integer nSQRSum, iSym, i, nData, iVer, j, N2, N1, iDSpin, nEle, iTmp, nEle1, nEle2, mTmp, iOff, n, iBas, iRC
real*8 GAP, eAlpha, eBeta, tmp

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
nSqrSum = 0
do iSym=1,nSym
  nSqrSum = nSqrSum+nBas(iSym)*nBas(iSym)
end do
!----------------------------------------------------------------------*
! Check start orbital priority list                                    *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(6,*) 'VecFind: LstVec',LstVec
#endif
do i=1,nStOpt
  Found = .false.
  if (LstVec(i) == 0) then
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Old scf orbitals'
#   endif
    call qpg_darray('SCF orbitals',Found,nData)
    if (Found .and. (nData == nSqrSum)) then
#     ifdef _DEBUGPRINT_
      write(6,*) 'Found orbitals'
#     endif
      call qpg_darray('OrbE',Found,nData)
      if (Found) then
#       ifdef _DEBUGPRINT_
        write(6,*) 'Found energies'
#       endif
        InVec = 8
        goto 101
      end if
    end if
  else if (LstVec(i) == 1) then
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Guessorb orbitals'
#   endif
    call qpg_darray('Guessorb',Found,nData)
    if (Found .and. (nData == nSqrSum)) then
      call qpg_darray('Guessorb energies',Found,nData)
      if (Found) then
        InVec = 9
        goto 101
      end if
    end if
  else if (LstVec(i) == 2) then
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Lumorb orbitals'
#   endif
    call F_Inquire(SCF_FileOrb,Found)
    if (Found) then
      if (isHDF5) then
#       ifdef _HDF5_
        iVer = -1
        call mh5_fetch_attr(fileorb_id,'NSYM',mynSym)
        call mh5_fetch_attr(fileorb_id,'NBAS',mynBas)
#       endif
      else
        call ChkVec(SCF_FileOrb,iVer,mynSym,mynBas,mynOrb,InfoLbl,iRc)
      end if
      if (iVer == 0) Found = .false.
      if (mynSym /= nSym) Found = .false.
      do j=1,nSym
        if (mynBas(j) /= nBas(j)) Found = .false.
      end do
    end if
    if (Found) then
#     ifdef _DEBUGPRINT_
      write(6,*) 'Found INPORB'
#     endif
      InVec = 2
      goto 101
    end if
  else if (LstVec(i) == 3) then
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Density'
#   endif
    InVec = 0
  else if (LstVec(i) == 4) then
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Core orbitals'
#   endif
    InVec = 0
    Found = .true.
    goto 101
  else
#   ifdef _DEBUGPRINT_
    write(6,'(2i5,a)') i,LstVec(i),' Die'
#   endif
    InVec = 0
    Found = .false.
    goto 101
  end if
# ifdef _DEBUGPRINT_
  write(6,*) 'LstVec(i),Found:',LstVec(i),Found
# endif
end do
Found = .false.
101 continue
#ifdef _DEBUGPRINT_
write(6,*) 'VecFind: InVec, Found=',InVec,Found
#endif
!----------------------------------------------------------------------*
! Did we find the requested orbitals?                                  *
!----------------------------------------------------------------------*
if (.not. Found) then
  cList = ''
  n2 = 0
  do i=1,nStOpt
    if (LstVec(i) == 0) then
      n1 = n2+1
      n2 = n1+9
      write(cList(n1:n2),'(a)') 'Old Scf, '
    else if (LstVec(i) == 1) then
      n1 = n2+1
      n2 = n1+10
      write(cList(n1:n2),'(a)') 'Guessorb, '
    else if (LstVec(i) == 2) then
      n1 = n2+1
      n2 = n1+8
      if (isHDF5) then
        write(cList(n1:n2),'(a)') 'HDF5, '
      else
        write(cList(n1:n2),'(a)') 'Lumorb, '
      end if
    else if (LstVec(i) == 3) then
      n1 = n2+1
      n2 = n1+13
      write(cList(n1:n2),'(a)') 'Old density, '
    else if (LstVec(i) == 4) then
      n1 = n2+1
      n2 = n1+6
      write(cList(n1:n2),'(a)') 'Core, '
    else if (LstVec(i) == 5) then
      n1 = n2+1
      n2 = n1+6
      write(cList(n1:n2),'(a)') 'NDDO, '
    else if (LstVec(i) == -1) then
      goto 290
    else
      n1 = n2+1
      n2 = n1+9
      write(cList(n1:n2),'(a)') 'Unknown, '
    end if
  end do
290 continue
  call SysAbendMsg('SCF:','Cannot find start orbitals according to list:',cList)
end if
!----------------------------------------------------------------------*
! What are the defaults for the different cases?                       *
!----------------------------------------------------------------------*
call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
Tot_El_Charge = Tot_Charge-Tot_Nuc_Charge
if (InVec == 0) then

  ! We will use core diagonalization

# ifdef _DEBUGPRINT_
  write(6,*) 'Using core diagonalization'
# endif
  OccSet = .false.
  FermSet = .true.
else if (Invec == 1) then

  ! We will use NDDO orbitals, should not be used!

# ifdef _DEBUGPRINT_
  write(6,*) 'Using NDDO orbitals'
# endif
else if (Invec == 2) then

  ! We will use Lumorb orbitals

# ifdef _DEBUGPRINT_
  write(6,*) 'Using Lumorb orbitals'
# endif
  call ChkLumo(OccSet,FermSet,SpinSet)
else if (Invec == 3) then

  ! We will use density as start, does it even work?

# ifdef _DEBUGPRINT_
  write(6,*) 'Using density'
# endif
else if (Invec == 4) then

  ! This is a restart case

# ifdef _DEBUGPRINT_
  write(6,*) 'Using Restart'
# endif
else if (Invec == 8) then

! We will use old SCF orbitals

# ifdef _DEBUGPRINT_
  write(6,*) 'Using SCF orbitals'
  write(6,*) 'tot_charge',tot_charge
  write(6,*) 'tot_el_charge',tot_el_charge
  write(6,*) 'tot_nuc_charge',tot_nuc_charge
# endif
  call qpg_iarray('SCF nOcc',Found,nData)
  idspin = 0
  if (Found) then
#   ifdef _DEBUGPRINT_
    write(6,*) 'vecfind: Alright, old scf orbitals it is'
#   endif
    nEle = 0
    if (nD == 1) then
      call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
      call Get_iScalar('SCF mode',iTmp)
      if (iTmp == 0) then
#       ifdef _DEBUGPRINT_
        write(6,*) 'Starting RHF with RHF orbitals'
#       endif
        do iSym=1,nSym
          nEle = nEle+2*nOcc(iSym,1)
        end do
      else
#       ifdef _DEBUGPRINT_
        write(6,*) 'Starting RHF with UHF orbitals'
#       endif
        call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
        nEle1 = 0
        nEle2 = 0
        do iSym=1,nSym
          nEle1 = nEle1+2*nOcc(iSym,1)
          nEle2 = nEle2+nOcc(iSym,1)+nOcc(iSym,2)
          idspin = idspin+nOcc(iSym,1)-nOcc(iSym,2)
        end do
        if (nEle1 /= nEle2) then
          !Tot_Charge = Tot_Nuc_Charge-nEle1
          !Tot_El_Charge = -nEle1
          CharSet = .true.
        end if
        nEle = nEle2
#       ifdef _DEBUGPRINT_
        write(6,*) 'After strange code'
        write(6,*) 'tot_charge',tot_charge
        write(6,*) 'tot_el_charge',tot_el_charge
        write(6,*) 'tot_nuc_charge',tot_nuc_charge
#       endif
      end if
    else
#     ifdef _DEBUGPRINT_
      write(6,*) 'Starting UHF with RHF/UHF orbitals'
#     endif
      call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
      call qpg_iarray('SCF nOcc_ab',Found,nData)
      if (Found) then
        call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
      else
        call Get_iarray('SCF nOcc',nOcc(1,2),nSym)
      end if
      do iSym=1,nSym
        nEle = nEle+nOcc(iSym,1)+nOcc(iSym,2)
        idspin = idspin+nOcc(iSym,1)-nOcc(iSym,2)
      end do
    end if
#   ifdef _DEBUGPRINT_
    write(6,*) 'idspin',idspin
    write(6,*) 'iAu_ab',iAu_ab
#   endif
    idspin = idspin-iAu_ab
    if ((abs(Tot_El_Charge+nEle) > Half) .or. (idspin /= 0)) then
      if (abs(Tot_El_Charge+nEle) > Half) then
#       ifdef _DEBUGPRINT_
        write(6,*) 'System have changed charge!'
#       endif
      end if
      if (idspin /= 0) then
#       ifdef _DEBUGPRINT_
        write(6,*) 'System have changed spin!'
#       endif
      end if
      if (CharSet .or. (idspin /= 0)) then
        OccSet = .false.
        FermSet = .true.
      else
        OccSet = .true.
        FermSet = .false.
        !Tot_Charge = Tot_Nuc_Charge-nEle
        !Tot_El_Charge = -nEle
      end if
    else
#     ifdef _DEBUGPRINT_
      write(6,*) 'System have same spin and charge'
#     endif
      OccSet = .true.
      FermSet = .false.
    end if
#   ifdef _DEBUGPRINT_
    write(6,*) 'OccSet  ',OccSet
    write(6,*) 'FermSet ',FermSet
    write(6,*) 'CharSet ',CharSet
    write(6,*) 'SpinSet ',SpinSet
    write(6,*) 'nOcc',nOcc
#   endif
  else
    OccSet = .false.
    FermSet = .true.
  end if
else if (Invec == 9) then

  ! We will use Guessorb orbitals

# ifdef _DEBUGPRINT_
  write(6,*) 'Using Guessorb orbitals'
# endif
  if (OccSet) then
#   ifdef _DEBUGPRINT_
    write(6,*) 'Occupation is set'
#   endif
  else if (FermSet) then
#   ifdef _DEBUGPRINT_
    write(6,*) 'Fermi is set'
#   endif
  else
#   ifdef _DEBUGPRINT_
    write(6,*) 'Must decide if to use Fermi'
#   endif
    if (nAufb(1) == -1) then
      mtmp = int(-Tot_El_Charge+0.1d0)
      if (nD == 1) then
        if (mod(mtmp,2) /= 0) then
          write(6,*) 'VecFind: Error in number of electrons'
          write(6,*) '         An even number of electrons ','are required by RHF, use UHF'
          write(6,*)
          call Abend()
        end if
        nAufb(1) = mtmp/2
      else
        nAufb(2) = (mtmp-iAu_ab)/2
        nAufb(1) = int(-Tot_El_Charge-nAufb(2))
      end if
    end if
#   ifdef _DEBUGPRINT_
    write(6,*) 'nAufb',nAufb
    write(6,*) 'Now figure out homo-lumo gap'
#   endif
    call qpg_darray('Guessorb energies',Found,nData)
    call mma_allocate(EOrb,nData,Label='EOrb')
    call get_darray('Guessorb energies',Eorb,nData)
    if (nD == 1) then
      call GetGap(Eorb,nData,nAufb(1),Gap,Ealpha)
    else
      call GetGap(Eorb,nData,nAufb(1),tmp,Ealpha)
      call GetGap(Eorb,nData,nAufb(2),Gap,Ebeta)
      Gap = min(tmp,Gap)
    end if
    call get_darray('Guessorb energies',Eorb,nData)
    if (Gap >= Half) then
      if (nD == 1) then
        iOff = 0
        do iSym=1,nSym
          n = 0
          do iBas=1,nBas(iSym)
            if (EOrb(iOff+iBas) < Ealpha) n = n+1
          end do
          nOcc(iSym,1) = n
          iOff = iOff+nBas(iSym)
        end do
      else
        iOff = 0
        do iSym=1,nSym
          n = 0
          do iBas=1,nBas(iSym)
            if (EOrb(iOff+iBas) < Ealpha) n = n+1
          end do
          nOcc(iSym,1) = n
          iOff = iOff+nBas(iSym)
        end do
        iOff = 0
        do iSym=1,nSym
          n = 0
          do iBas=1,nBas(iSym)
            if (EOrb(iOff+iBas) < Ebeta) n = n+1
          end do
          nOcc(iSym,2) = n
          iOff = iOff+nBas(iSym)
        end do
      end if
      OccSet = .true.
      FermSet = .false.
#     ifdef _DEBUGPRINT_
      write(6,*) 'Decided on occupation list'
#     endif
    else
      OccSet = .false.
      FermSet = .true.
#     ifdef _DEBUGPRINT_
      write(6,*) 'Decided on Fermi aufbau'
#     endif
    end if
    call mma_deallocate(EOrb)
#   ifdef _DEBUGPRINT_
    write(6,*) 'Gap is',Gap
#   endif
  end if
else

  ! This case should not appear

  call SysAbendMsg('SCF:','Internal error in VecFind!','InVec have illegal value')
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(6,'(a,i2)') 'VecFind: InVec=',InVec
write(6,*) 'OccSet=',OccSet
write(6,*) 'FermSet=',FermSet
#endif
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine VecFind
