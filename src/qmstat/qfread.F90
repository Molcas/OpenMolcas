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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  Qfread
!
!> @brief
!>   Read in all data that comes from external Molcas routines and prepare various quantities, such as MME
!> @author A. Ohrn
!>
!> @details
!> This subroutine handles the interaction with the rest of Molcas.
!> Here orbitals and various matrices are stored and to some extent
!> modified to our purpose. We call on this subroutine even when
!> we are running classical stuff. The reason for this is that we
!> wish to collect some numbers, but really, we could easily have
!> constructed thing differently so that this subroutine only would
!> be called when quantum-classical stuff is running. For RASSI
!> implementation, we also read in and transform the transition density
!> matrix. At the end, we also call the routines that make the
!> multicenter multipole expansion.
!>
!> @note
!> Seward is mandatory for both SCF and RASSI. For SCF also,
!> Motra, Averd; for RASSI also, RASSCF and RASSI.
!>
!> @param[in]  iQ_Atoms
!> @param[out] nAtomsCC Atoms in solvent
!> @param[out] Coord    Unique coordinates of the atoms in the molecule in the QM-region
!> @param[out] nBas     Number of basis functions in QM-region
!> @param[out] nBasCC   Like nBas but for a solvent molecule
!> @param[out] nOcc     The total number of basis functions that belong to a certain basis-function type.
!> @param[out] natyp    Number of atoms of the i:th basis-function type
!> @param[out] nntyp    Number of basis-function types.
!***********************************************************************

!******JoseMEP the last three variables are included to the MEP calculation
subroutine Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBasCC,nOcc,natyp,nntyp)

use qmstat_global, only: Alfa, ATitle, BasOri, Beta, c_orbene, ChaNuc, Cont, CT, Dont, info_atom, iOrb, iPrint, iQang, iQn, &
                         iWoGehenC, iWoGehenQ, Joblab, lmax, mPrimus, MxAngqNr, MxSymQ, nBA_C, nBA_Q, nBonA_C, nBonA_Q, nCBoA_C, &
                         nCBoA_Q, nCnC_C, nPrimus, PotNuc, QmType, SavOri, Trans, V1, V3
use qmstat_procedures, only: GiveMeInfo
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms
integer(kind=iwp), intent(out) :: nAtomsCC, nBas(MxSymQ), nBasCC(1), nOcc(iQ_Atoms), natyp(iQ_Atoms), nntyp
real(kind=wp), intent(out) :: Coord(3,iQ_Atoms)
integer(kind=iwp) :: i, iAtom, iBas, icont, iDummy(1), iErr, iLu, ind, indold, iold, iWarn, ix, j, jnd, k, kk, kold, l, lLine, na, &
                     nnaa, nntypC, nSize, nSym, nSymCC, ntBas
real(kind=wp) :: ChgeCC(3), CoordCC(3*3), Dummy(1)
character(len=120) :: Line, StLine
character(len=100) :: OrbName, Title
character(len=10) :: WhatGet
integer(kind=iwp), allocatable :: iC_Icon(:,:), Icon(:,:), natypC(:), nfSh(:,:), nSh(:)
real(kind=wp), allocatable :: C(:,:), Chge(:), Cmo(:,:), Cmo_S(:), E(:,:), Occu(:), Oe(:), Tmp(:,:), TransC(:)
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

!----------------------------------------------------------------------*
! Enter                                                                *
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
! Print the joblabel. It is obtained in get_input                      *
!----------------------------------------------------------------------*
if (ATitle) then
  lLine = len(Line)
  do i=1,lLine
    StLine(i:i) = '*'
  end do
  write(u6,*)
  do i=1,6
    Line = ''
    if ((i == 1) .or. (i == 6)) Line = StLine
    if (i == 3) Line = 'Project:'
    if (i == 4) write(Line,'(A72)') JobLab
    call Center_Text(Line)
    write(u6,*) '*'//Line//'*'
  end do
  write(u6,*)
end if
write(u6,*) 'Auxiliary data being read and pre-processed.'
!----------------------------------------------------------------------*
! Collect some data from RUNFILE about the QM-region molecule.         *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
if (nSym /= 1) then !A minor restriction, no symmetry allowed, i.e. nSym=1.
  write(u6,*)
  write(u6,*) ' QmStat does not run with symmetry!'
  write(u6,*) ' The perturbation from the solvent breaks all symmetry.'
  call Quit(_RC_GENERAL_ERROR_)
end if
!call Get_iScalar('Unique atoms',iQ_Atoms)
call mma_allocate(Chge,iQ_Atoms,label='Chge')
call Get_dArray('Nuclear charge',Chge,iQ_Atoms)
call Get_dArray('Unique Coordinates',Coord,3*iQ_Atoms)
call Get_dArray('Center of Mass',CT,3)
call Get_iArray('nBas',nBas,nSym)

!----------------------------------------------------------------------*
! Print elementary information about molecule.                         *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,*) '     ------------------------------'
write(u6,*) '     |       QM-region data       |'
write(u6,*) '     ------------------------------'
write(u6,*)
write(u6,'(A,15X,I5)') '      Number of basis functions:',(nBas(i),i=1,nSym)
write(u6,'(A,3F10.6)') '      Centre of mass =',(CT(kk),kk=1,3)
call PrCoor()

!----------------------------------------------------------------------*
! If the Qmtype is SCF we now want the orbitals. On the other hand if  *
! we are running RASSI, another route must be taken, hence here we     *
! inquire which QM-method is used.                                     *
!----------------------------------------------------------------------*
if (QmType(1:3) == 'SCF') then
  !--------------------------------------------------------------------*
  !    SSS    CCC  FFFFFF                                              *
  !   SS    CC     FF                                                  *
  !    SS  CC      FFFF                                                *
  !     SS  CC     FF                                                  *
  !   SSS    CCC   FF                                                  *
  !--------------------------------------------------------------------*

  !--------------------------------------------------------------------*
  ! Print information about orbitals and store coefficients in new     *
  ! variable.                                                          *
  !--------------------------------------------------------------------*
  iLu = 15
  iLu = IsFreeUnit(iLu)
  write(OrbName,'(A)') 'AVEORB'
  write(WhatGet,'(A)') 'CO'
  iWarn = 1
  call mma_allocate(Cmo,nBas(1),nBas(1),label='Cmo')
  call mma_allocate(Occu,nBas(1),label='Occu')
  call RdVec(OrbName,iLu,WhatGet,nSym,nBas,nBas,Cmo,Occu,Dummy,iDummy,Title,iWarn,iErr)
  if (iErr /= 0) then
    write(u6,*)
    write(u6,*) 'Error when reading AVEORB'
    write(u6,*)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  call mma_allocate(V1,nBas(1),iOrb(1),label='OrbCoeffQ')
  V1(:,:) = Cmo(:,1:iOrb(1))
  !write(u6,'(A,I4)') '      Number of Orbitals:',iOrb(1)
  call Get_dScalar('PotNuc',PotNuc)

else if (QmType(1:4) == 'RASS') then
  !--------------------------------------------------------------------*
  !   RRRR      AA      SSS   SSS   II                                 *
  !   RR  R    A  A    SS    SS     II                                 *
  !   RR R    AAAAAA    SS    SS    II                                 *
  !   RRR     AA  AA     SS    SS   II                                 *
  !   RR R    AA  AA   SSS   SSS    II                                 *
  !--------------------------------------------------------------------*
  call TdmTrans(nBas)
end if
write(u6,*)
write(u6,*)
!----------------------------------------------------------------------*
! Compute various information about system. This we use for computing  *
! integrals later. GiveMeInfo collects stuff from seward, sometime with*
! some recomputations.                                                 *
!----------------------------------------------------------------------*
call mma_allocate(nSh,iQ_Atoms,label='nSh')
call mma_allocate(nfSh,iQ_Atoms,MxAngqNr,label='nfSh')
call mma_allocate(nBA_Q,iQ_Atoms,label='nBA_Q')
call mma_allocate(nBonA_Q,iQ_Atoms,label='nBonA_Q')
call mma_allocate(nCBoA_Q,iQ_Atoms,MxAngqNr,label='nCBoA_Q')
call GiveMeInfo(nntyp,natyp,BasOri,Icon,nPrimus,nBA_Q,nCBoA_Q,nBonA_Q,E,C,nSh,nfSh,nSize,iPrint,iQ_Atoms,MxAngqNr,Trans,ntBas)
iBas = 0
iAtom = 0
kold = 1
iold = 1
indold = 0
nOcc(1:nntyp) = 0
call mma_allocate(alfa,ntBas,0,label='alfa')
call mma_allocate(cont,ntBas,0,label='cont')
call mma_allocate(iWoGehenQ,ntBas,2*MxAngqNr-1,label='iWoGehenQ')
call mma_allocate(iQang,ntBas,label='iQang')
do i=1,nntyp
  na = natyp(i)
  do j=1,na
    ind = 0
    jnd = 0
    iAtom = iAtom+1
    ChaNuc(iAtom) = Chge(iAtom)
    info_atom(iAtom) = int(Chge(iAtom))
    do k=1,nSh(i)
      nnaa = nfSh(i,k)
      do l=1,nnaa
        ibas = ibas+1
        indold = indold+1
        nOcc(i) = nOcc(i)+2*k-1
        ind = ind+1
        iQang(ibas) = k
        icont = Icon(i,ind)
        !iCharOnBasQ(ibas) = int(Chge(iAtom))
        ! Here we construct an array of indices which is used to put right AO-overlap in right matrix pos.
        do ix=1,2*k-1
          if (k /= kold) then
            if (i /= iold) then
              Indold = Indold+nfSh(iold,kold)*(2*kold-2)
              iold = i
            else
              Indold = Indold+nfSh(i,kold)*(2*kold-2)
            end if
            kold = k
          end if
          iWoGehenQ(ibas,ix) = indold+nnaa*(ix-1)
        end do
        if (icont > size(cont,2)) then
          call mma_allocate(Tmp,size(alfa,1),icont,label='Tmp')
          Tmp(:,1:size(alfa,2)) = alfa
          call mma_deallocate(alfa)
          call move_alloc(Tmp,alfa)
          call mma_allocate(Tmp,size(cont,1),icont,label='Tmp')
          Tmp(:,1:size(cont,2)) = cont
          call mma_deallocate(cont)
          call move_alloc(Tmp,cont)
        end if
        alfa(ibas,1:icont) = E(i,jnd+1:jnd+icont)
        cont(ibas,1:icont) = C(i,jnd+1:jnd+icont)
        jnd = jnd+icont
      end do
    end do
  end do
end do
! Now we do not need them, so deallocate
call mma_deallocate(Chge)
call mma_deallocate(Icon)
call mma_deallocate(E)
call mma_deallocate(C)

!----------------------------------------------------------------------*
! Obtain and print information about solvent. This requires a renaming *
! of the runfile.                                                      *
!----------------------------------------------------------------------*
call NameRun('WRUNFIL')
call Get_iScalar('nSym',nSymCC)
if (nSymCC /= 1) then
  write(u6,*)
  write(u6,*) ' QmStat does not run with symmetry!'
  call Quit(_RC_GENERAL_ERROR_)
end if
call Get_iScalar('Unique atoms',nAtomsCC)
if (nAtomsCC /= 3) then
  write(u6,*)
  write(u6,*) 'Now now... what strange solvent molecule do you try to feed QmStat with?'
  call Quit(_RC_GENERAL_ERROR_)
end if
call Get_dArray('Nuclear charge',ChgeCC,nAtomsCC)
call Get_dArray('Unique Coordinates',CoordCC,3*nAtomsCC)
call Get_iArray('nBas',nBasCC,nSymCC)
write(u6,*)
write(u6,*) '     ------------------------------'
write(u6,*) '     |   Solvent molecule data     |'
write(u6,*) '     ------------------------------'
write(u6,*)
write(u6,'(A,15X,I5)') '      Number of basis functions:',(nBasCC(i),i=1,nSymCC)
call PrCoor()

! Collect information about the solvent orbitals.

iLu = IsFreeUnit(16)
write(OrbName,'(A)') 'SOLORB'
write(WhatGet,'(A)') 'CE'
iWarn = 1
call mma_allocate(Oe,sum(nBasCC),label='OrbitalEnergy')
call mma_allocate(Cmo_S,nBasCC(1)**2,label='Cmo_S')
call RdVec(OrbName,iLu,WhatGet,nSymCC,nBasCC,nBasCC,Cmo_S,Dummy,Oe,iDummy,Title,iWarn,iErr)
call mma_allocate(c_orbene,iOrb(2),label='c_orbene')
c_orbene(1:iOrb(2)) = Oe(1:iOrb(2))
call mma_deallocate(Oe)

! We should not need two solvent orbital vectors, so this should
! be removed when the orbital rotation routine is fixed.

call mma_allocate(V3,nBasCC(1),iOrb(2),label='V3')
do j=1,iOrb(2)
  do k=1,nBasCC(1)
    V3(k,j) = Cmo_S(k+(j-1)*nBasCC(1))
  end do
end do
call mma_deallocate(Cmo_S)
!write(u6,'(A,I4)') '      Number of Orbitals:',iOrb(2)
write(u6,*)
write(u6,*)
!----------------------------------------------------------------------*
! And now basis set information.                                       *
!----------------------------------------------------------------------*
if (nAtomsCC /= iQ_Atoms) then
  call mma_deallocate(nSh)
  call mma_deallocate(nfSh)
  call mma_allocate(nSh,nAtomsCC,label='nSh')
  call mma_allocate(nfSh,nAtomsCC,MxAngqNr,label='nfSh')
end if
call mma_allocate(natypC,nAtomsCC,label='natypC')
call mma_allocate(nBA_C,nAtomsCC,label='nBA_C')
call mma_allocate(nBonA_C,nAtomsCC,label='nBonA_C')
call mma_allocate(nCBoA_C,nAtomsCC,MxAngqNr,label='nCBoA_C')
call GiveMeInfo(nntypC,natypC,SavOri,iC_Icon,mPrimus,nBA_C,nCBoA_C,nBonA_C,E,C,nSh,nfSh,nSize,iPrint,nAtomsCC,MxAngqNr,TransC,ntBas)
iBas = 0
iAtom = 0
kold = 1
iold = 1
indold = 0
call mma_allocate(beta,ntBas,0,label='beta')
call mma_allocate(dont,ntBas,0,label='dont')
call mma_allocate(iWoGehenC,ntBas,2*MxAngqNr-1,label='iWoGehenC')
call mma_allocate(iQn,ntBas,label='iQn')
call mma_allocate(nCnC_C,ntBas,label='nCnC_C')
Lmax = ntBas
do i=1,nntypC !Like the corresponding thing above for the QM-region.
  na = natypC(i)
  do j=1,na
    ind = 0
    jnd = 0
    iAtom = iAtom+1
    do k=1,nSh(i)
      nnaa = nfSh(i,k)
      do l=1,nnaa
        ibas = ibas+1
        indold = indold+1
        nCnC_C(ibas) = nnaa
        ind = ind+1
        iQn(ibas) = k
        icont = iC_Icon(i,ind)
        !iCharOnBasC(ibas) = int(ChgeCC(iAtom))
        ! Here we construct an array of indices which is used to put right AO-overlap in right matrix pos.
        do ix=1,2*k-1
          if (k /= kold) then
            if (i /= iold) then
              Indold = Indold+nfSh(iold,kold)*(2*kold-2)
              iold = i
            else
              Indold = Indold+nfSh(i,kold)*(2*kold-2)
            end if
            kold = k
          end if
          iWoGehenC(ibas,ix) = indold+nnaa*(ix-1)
        end do
        if (icont > size(dont,2)) then
          call mma_allocate(Tmp,size(beta,1),icont,label='Tmp')
          Tmp(:,1:size(beta,2)) = beta
          call mma_deallocate(beta)
          call move_alloc(Tmp,beta)
          call mma_allocate(Tmp,size(dont,1),icont,label='Tmp')
          Tmp(:,1:size(dont,2)) = dont
          call mma_deallocate(dont)
          call move_alloc(Tmp,dont)
        end if
        beta(ibas,1:icont) = E(i,jnd+1:jnd+icont)
        dont(ibas,1:icont) = C(i,jnd+1:jnd+icont)
        jnd = jnd+icont
      end do
    end do
  end do
end do
if (size(TransC) > size(Trans)) then
  call mma_deallocate(Trans)
  call move_alloc(TransC,Trans)
else
  call mma_deallocate(TransC)
end if
! Now we do not need them, so deallocate.
call mma_deallocate(natypC)
call mma_deallocate(nSh)
call mma_deallocate(nfSh)
call mma_deallocate(iC_Icon)
call mma_deallocate(E)
call mma_deallocate(C)

!----------------------------------------------------------------------*
! The multipoles and the Hamiltonian matrix are radically different    *
! between the QM-method alternatives, so once more an inquire.         *
!----------------------------------------------------------------------*
call NameRun('#Pop')
if (QmType(1:3) == 'SCF') then
  call ScfHandM(Cmo,nBas,iQ_Atoms,nOcc,natyp,nntyp,Occu)
  call mma_deallocate(Cmo)
  call mma_deallocate(Occu)
else if (QmType(1:4) == 'RASS') then
  call RassiHandM(nBas,iQ_Atoms,nOcc,natyp,nntyp)
end if

!----------------------------------------------------------------------*
! Here is the end.                                                     *
!----------------------------------------------------------------------*
return

end subroutine Qfread
