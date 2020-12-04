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
! Copyright (C) 2009, Giovanni Li Manni                                *
!               2009, Francesco Aquilante                              *
!               2020, Oskar Weser                                      *
!***********************************************************************

module desymmetrize_mod
    use Basis_Info, only: nBas, nCnttp, dbsc, Shells, MolWgh
    use Center_Info, only: dc
    use definitions, only: wp
    use linalg_mod, only: abort_, verify_
    use Symmetry_Info, only: nIrrep, lIrrep
    use stdalloc, only: mma_allocate, mma_deallocate
    use sorting, only: swap, sort, argsort
    use sorting_funcs, only: leq_r, geq_r
    use info_expbas_mod

    implicit none
    private
    public :: desym

    integer, parameter :: nNumber = 61
    character(len=1), parameter :: &
        number(nNumber) = &
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', &
         'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', &
         'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', &
         'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', &
         'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', &
         'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', &
         'Z']
    !> Different kinds of orbitals
    !> f, i, 1, 2, 3, s, d
    integer, parameter :: n_orb_kinds = 7


! NOTE: These global variables are ugly as hell, but we need
!  it to support the shitty SUN compiler.
#ifndef INTERNAL_PROC_ARG
    !> This is the orbital kind for each orbital.
    integer, allocatable :: kind_per_orb(:)
    real(wp), allocatable :: energy(:), occ(:)
    real(wp), allocatable :: CMO(:, :)
#endif

contains

! symmetry-----> C1 INPORB
    Subroutine desym()
!***********************************************************************
!                                                                      *
! Purpose: Convert INPORB file with symmetry to DeSymOrb               *
!          file with C1 symmetry.                                      *
!                                                                      *
!          G. Li Manni, F. Aquilante, University of Geneva, June 2009. *
!          Oskar Weser, MPI FKF Stuttgart,  November 2020.             *
!                                                                      *
! Example of input:                                                    *
!                                                                      *
!   &SEWARD                                                            *
!    coord                                                             *
!     file.xyz                                                         *
!    basis                                                             *
!     ano-rcc-mb                                                       *
!    oneonly                                                           *
!    nodk                                                              *
!                                                                      *
!   >>COPY $CurrDir/scf_sym.ScfOrb INPORB                              *
!                                                                      *
!   &EXPBAS                                                            *
!    NoEx                                                              *
!    Desy                                                              *
!   END                                                                *
!                                                                      *
!***********************************************************************
#include "Molcas.fh"
#include "WrkSpc.fh"
        real(wp), parameter :: EorbThr = 50._wp
        real(wp) :: Coor(3, MxAtom), Znuc(MxAtom)
        character(len=LENIN) :: AtomLabel(MxAtom)
        character(len=512) :: FilesOrb
        character(len=LENIN8), allocatable :: label(:)
        character(len=8) :: MO_Label(maxbfn)
        integer :: ibas_lab(MxAtom), nOrb(8)
#ifdef INTERNAL_PROC_ARG
        !> This is the orbital kind for each orbital.
        integer, allocatable :: kind_per_orb(:)
        real(wp), allocatable :: energy(:), occ(:)
        real(wp), allocatable :: CMO(:, :)
#endif
        !> This is the number of orbitals for every kind.
        integer :: n_kinds(n_orb_kinds)

        character(len=LENIN8 + 1) :: gtolabel(maxbfn)
        character(len=50) :: VTitle
        character(len=128) :: SymOrbName
        logical :: exists, found
        logical, parameter :: y_cart = .false.


        integer :: nAtom, nData, nDeg, nTot, nTot2, nB
        integer :: iCnttp, iAngMx_Valence
        integer :: ierr
        integer :: ipCent, ipCent2, ipCent3
        integer :: ipPhase, ipC2, ipV
        integer :: mAdOcc, mAdEor, mAdCMO
        integer :: file_id, iWfType
        integer, parameter :: notSymm = 1, arbitrary_number = 42, noUHF = 0, iWarn = 1
        integer :: iatom, iDeg, ishell, iIrrep

        integer :: mdc, kk, i, j, ik, k, l, kk_Max, ii, iB, ipp, ic, iv
        integer :: ipc, icontr, nBasisi, icntr

        integer, save :: iRc = 0

        file_id = arbitrary_number


        Call f_Inquire('RUNFILE', exists)
        call verify_(exists, 'Error finding RUNFILE')

        !-----Read the characteristics of all different basis sets,
        !     provide each atom with a nuclear charge and establish
        !     a link between an atom and its basis set ---
        !
        !     NOTICE!!!
        !     This call will also fill info.fh and the dynamic storage in
        !     Work(ipInf)
        !
        Call Inter1(AtomLabel, iBas_Lab, Coor, Znuc, nAtom)
        Call Qpg_iArray('nOrb', found, nData)
        If (found) Then
            Call Get_iArray('nOrb', nOrb, nData)
        Else
            nOrb(:nIrrep) = nBas(:nIrrep)
        End If

        iAngMx_Valence = maxval(dbsc%nVal - 1, mask=.not. (dbsc%Aux .or. dbsc%Frag))

        !     Compute memory requirements and allocate memory
        !
        nB = sum(nBas(0:nIrrep - 1))
        Call GetMem('ICENT', 'ALLO', 'INTE', ipCent, 8 * nB)
        Call GetMem('IPHASE', 'ALLO', 'INTE', ipPhase, 8 * nB)
        Call GetMem('nCENT', 'ALLO', 'INTE', ipCent2, nB)
        Call GetMem('ICENTER', 'ALLO', 'INTE', ipCent3, nB)
        Call GetMem('CMO2', 'ALLO', 'REAL', ipC2, nB**2)
        Call GetMem('VECTOR', 'ALLO', 'REAL', ipV, nB**2)
        call dcopy_(nB**2, [0._wp], 0, Work(ipV), 1)
        Call FZero(Work(ipC2), nB**2)

        !     Read exponents and contraction coefficients of each unique basis.
        !     Write the present basis set (iCnttp) to the molden.input file for
        !     the appropriate atoms.
        !     Moreover, a list is constructed which contains a label for each
        !     GTO (gtolabel). This list follows the MOLDEN order of GTOs.
        !     Later this list will be used in the transformation of sabf (the
        !     symmetry adapted basis functions).
        !
        iatom = 0; mdc = 0; kk = 0
        Do iCnttp = 1, nCnttp             ! loop over unique basis sets
            If (.not. (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) then
                Do iCntr = 1, dbsc(iCnttp)%nCntr! loop over symmetry unique centers
                    mdc = mdc + 1
                    nDeg = nIrrep / dc(mdc)%nStab
                    Do iDeg = 1, nDeg             ! loop over centers
                        iAtom = iAtom + 1
                        !
                        If (dbsc(iCnttp)%nVal > 6) Then
                            Write (6, *) 'Desym: too high angular momentum!'
                            write (6, *) 'iCnttp and dbsc(iCnttp)%nVal= ' &
                                , iCnttp, dbsc(iCnttp)%nVal
                            Call Abend()
                        End If
                        !
                        Do l = 0, dbsc(iCnttp)%nVal - 1
                            ishell = dbsc(iCnttp)%iVal + l
                            nBasisi = Shells(iShell)%nBasis
                            !              write(6,*) 'nBasisi', Shells(iShell)%nBasis
                            If (nBasisi > nNumber) Then
                                Write (6, *) 'Desym: too many contracted functions!'
                                Write (6, *) 'nBasisi=', Shells(iShell)%nBasis
                                Call Abend()
                            End If

                            !             Iterate over each contracted GTO
                            !
                            If (l == 0) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'01s     '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If (l == 1) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02px    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02py    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02pz    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 2) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 3) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 4) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            EndIf
                            If ((l == 5) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h04-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h04+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                        End Do
                    End Do
                End Do
            end if
        End Do
        kk_Max = kk
        If (nB > kk_max) Then
            Write (6, *) 'nB > kk_max'
            Write (6, *) 'nB,kk_Max=', nB, kk_Max
            call ClsSew()
            return
        End If

        nTot = sum(nBas(0 : nIrrep - 1))
        nTot2 = sum(nBas(0 : nIrrep - 1)**2)

        call mma_allocate(kind_per_orb, nTot)
        Call GetMem('Occ', 'Allo', 'Real', mAdOcc, nTot)
        Call GetMem('Eor', 'Allo', 'Real', mAdEor, nTot)
        Call GetMem('CMO', 'Allo', 'Real', mAdCMO, nTot2)
        Call FZero(Work(mAdOcc), nTot)
        Call FZero(Work(mAdEor), nTot)
        Call FZero(Work(mAdCMO), nTot2)

        !---- Read HF CMOs from file
        FilesOrb = EB_FileOrb
        If (len_trim(FilesOrb) == 0) FilesOrb = 'INPORB'
        if (DoExpbas) FilesOrb = 'EXPORB'
        Call RdVec_(trim(FilesOrb), file_id, 'COEI', noUHF, nIrrep, nBas, nBas, &
                    Work(mAdCMO), Work(ip_Dummy), &
                    Work(mAdOcc), Work(ip_Dummy), &
                    Work(mAdEor), Work(ip_Dummy), &
                    kind_per_orb, VTitle, iWarn, iErr, iWfType)

        if (iErr /= 0) then
            iRc = 1
            return
        end if

        !     Get the coeff. of sym adapted basis functions (ipC2)
        !
        Call Dens_IF_SCF(Work(ipC2), Work(mAdCMO), 'F')
        Call GetMem('CMO', 'Free', 'Real', mAdCMO, nTot2)

        !      Back 'transformation' of the symmetry adapted basis functions.
        !      Probably somewhat clumsy, but it seems to work.If someone
        !      knows a more elegant way to do it, please improve this part!
        !
        !      PART 1: Obtain symmetry information (soout), construct a label
        !              for each sabf, which will be used in part 2 to find the
        !              corresponding GTO in the MOLDEN list by comparing with
        !              gtolabel
        !
        !      nB       --- Total number of contracted basis functions
        !      ipcent2  --- degeneracy of a basis function
        !      ipCent   --- centres over which the basis function is
        !                   delocalized
        !      ipPhase  --- phase of the AO in the linear combination
        !
        Call mma_allocate(Label, MaxBfn + MaxBfn_Aux, label='Label')
        Call icopy(8 * nB, [0], 0, iWork(ipPhase), 1)
        Call icopy(8 * nB, [0], 0, iWork(ipCent), 1)
        Call SOout(label, iWork(ipCent), iWork(ipPhase))
        ipc = 0
        Do iContr = 1, nB
            iWork(ipCent2 + iContr - 1) = 0
            Do k = 1, 8
                If (iWork(ipCent + ipc) /= 0) then
                    iWork(ipcent2 + iContr - 1) = iWork(ipCent2 + iContr - 1) + 1
                end if
                ipc = ipc + 1
            End Do
        End Do
        ! Part 2: -Take a MOLCAS symmetry functions (loop i)
        !         -Find the corresponding label in the MOLDEN list (loop j)
        !         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
        !          by the appropriate factor (ipPhase),and divide by the number of
        !          centres over which the sabf is delocalized (ipCent3).
        !         -The vectors are copied by rows!
        !
        !     loop over MOLCAS symmetry functions
        i = 0
        ik = 0
        Do iIrrep = 0, nIrrep - 1
            !
            Do iB = 1, nBas(iIrrep)
                i = i + 1
                !
                If (iB == 1) Then
                    ik = 1
                else
                    If (label(i - 1) == label(i)) Then
                        ik = ik + 1
                    else
                        ik = 1
                    End If
                End If
                Write (MO_Label(i), '(I5,A3)') iB, lirrep(iIrrep)
                !
                Do j = 1, nB
                    !
                    !
                    If (gtolabel(j) == label(i)//number(ik)) Then
                        Do k = 1, 8
                            ipc = (i - 1) * 8 + k - 1
                            ipp = ipc
                            If (iWork(ipCent + ipc) == iWork(ipcent3 + j - 1)) Then
                                Do ii = 1, nB
                                    ic = (ii - 1) * nB + (i - 1)
                                    iv = (ii - 1) * nB + (j - 1)
                                    If (MolWgh == 0) Then
                                        Work(ipV + iv) = Work(ipV + iv) &
                                                         + Work(ipC2 + ic) &
                                                         * DBLE(iWork(ipPhase + ipp)) &
                                                         / DBLE(iWork(ipcent2 + i - 1))
                                    Else
                                        Work(ipV + iv) = Work(ipV + iv) &
                                                         + Work(ipC2 + ic) &
                                                         * DBLE(iWork(ipPhase + ipp)) &
                                                         / Sqrt(DBLE(iWork(ipcent2 + i - 1)))
                                    End If
                                End Do
                            End If
                        End Do
                    End If
                End Do
            End Do
        End Do
        Call mma_deallocate(Label)



        !**************************** START SORTING ****************************

        call mma_allocate(CMO, nTot, nTot)
        call mma_allocate(occ, nTot)
        call mma_allocate(energy, nTot)

        energy(:) = Work(mAdEor : mAdEor + nTot - 1)
        occ(:) = Work(mAdOcc : mAdocc + nTot - 1)

        do i = 0, nTot - 1
            l = ipV + nTot * i
            CMO(:, i + 1) = work(l : l + nTot - 1)
        end do

#ifdef INTERNAL_PROC_ARG
        call reorder_orbitals(nTot, kind_per_orb, CMO, occ, energy)
#else
        call reorder_orbitals()
#endif

        n_kinds(:) = 0
        do i = lbound(kind_per_orb, 1), ubound(kind_per_orb, 1)
            n_kinds(kind_per_orb(i)) = n_kinds(kind_per_orb(i)) + 1
        end do


        SymOrbName = 'DESORB'
        VTitle = 'Basis set desymmetrized orbital file DESORB'
        Call WrVec_(SymOrbName, file_id, 'COEI', noUHF, notSymm, &
                    [nTot], [nTot], &
                    CMO, Work(ip_Dummy), &
                    occ, Work(ip_Dummy), &
                    energy, Work(ip_Dummy), &
                    n_kinds, VTitle, iWFtype)
        call Add_Info('desym CMO', CMO, 999, 8)

        call mma_deallocate(occ)
        call mma_deallocate(CMO)
        call mma_deallocate(energy)
        call mma_deallocate(kind_per_orb)
        Call GetMem('Eor', 'Free', 'Real', mAdEor, nTot)
        Call GetMem('Occ', 'Free', 'Real', mAdOcc, nTot)
        Call GetMem('ICENT', 'FREE', 'INTE', ipCent, 8 * nB)
        Call GetMem('IPHASE', 'FREE', 'INTE', ipPhase, 8 * nB)
        Call GetMem('nCENT', 'FREE', 'INTE', ipCent2, nB)
        Call GetMem('ICENTER', 'FREE', 'INTE', ipCent3, nB)
        Call GetMem('CMO2', 'FREE', 'REAL', ipC2, nB**2)
        Call GetMem('VECTOR', 'FREE', 'REAL', ipV, nB**2)

        call ClsSew()

    end subroutine

#ifdef INTERNAL_PROC_ARG

    subroutine reorder_orbitals(nTot, kind_per_orb, CMO, occ, energy)
        integer, intent(in) :: nTot
        integer, intent(inout) :: kind_per_orb(nTot)
        real(wp), intent(inout) :: CMO(nTot, nTot), occ(nTot), energy(nTot)

        integer :: i
        integer, allocatable :: idx(:)

        allocate(idx(nTot))

        idx(:) = [(i, i = 1, nTot)]

        call sort(idx, closure_compare)

        kind_per_orb(:) = kind_per_orb(idx)
        energy(:) = energy(idx)
        CMO(:, :) = CMO(:, idx)
        occ(:) = occ(idx)

        contains

            !> @brief
            !>  Comparison function for generic sort.
            !>
            !> @details
            !>  Sort non-strict i.e. (compare(i, j) .and. compare(j, i)) can be true)
            !>  Sort first by orbital kind ascendingly (frozen, inactive, RAS1, ...),
            !>      second by occupation number descendingly (2.0, 2.0, 1.x, 0., ...),
            !>      third by energy ascendingly (-3., -2., -2., 0., 1., ...).
            !>  Note, that `sort` uses a stable sorting algorithm and since the
            !>  input orbitals are automatically by irrep the output
            !>  will be sorted last ascendingly by irrep.
            pure function closure_compare(i, j) result(res)
                integer, intent(in) :: i, j
                logical :: res

                if (kind_per_orb(i) /= kind_per_orb(j)) then
                    res = kind_per_orb(i) < kind_per_orb(j)
                else if (occ(i) /= occ(j)) then
                    res = occ(i) > occ(j)
                else if (energy(i) /= energy(j)) then
                    res = energy(i) < energy(j)
                else
                    ! All relevant values are equal and our comparison has
                    ! to be non-strict.
                    res = .true.
                end if
            end function
    end subroutine

#else

    subroutine reorder_orbitals()
        integer :: nTot, i
        integer, allocatable :: idx(:)

        nTot = size(kind_per_orb)
        allocate(idx(nTot))
        idx(:) = [(i, i = 1, nTot)]

        call sort(idx, compare)

        kind_per_orb(:) = kind_per_orb(idx)
        energy(:) = energy(idx)
        CMO(:, :) = CMO(:, idx)
        occ(:) = occ(idx)

! There was a wrong and stupid warning from GFortran 4.8
! If we stop supporting this compiler remove this ugly workaround.
#ifdef _WARNING_WORKAROUND_
        if (.false.) then
            if (compare(1, 2)) continue
        end if
#endif
    end subroutine

    !> @brief
    !>  Comparison function for generic sort.
    !>
    !> @details
    !>  Sort non-strict i.e. (compare(i, j) .and. compare(j, i)) can be true)
    !>  Sort first by orbital kind ascendingly (frozen, inactive, RAS1, ...),
    !>      second by occupation number descendingly (2.0, 2.0, 1.x, 0., ...),
    !>      third by energy ascendingly (-3., -2., -2., 0., 1., ...).
    !>  Note, that `sort` uses a stable sorting algorithm and since the
    !>  input orbitals are automatically by irrep the output
    !>  will be sorted last ascendingly by irrep.
    pure function compare(i, j) result(res)
        integer, intent(in) :: i, j
        logical :: res

        if (kind_per_orb(i) /= kind_per_orb(j)) then
            res = kind_per_orb(i) < kind_per_orb(j)
        else if (occ(i) /= occ(j)) then
            res = occ(i) > occ(j)
        else if (energy(i) /= energy(j)) then
            res = energy(i) < energy(j)
        else
            ! All relevant values are equal and our comparison has
            ! to be non-strict.
            res = .true.
        end if
    end function

#endif

end module
