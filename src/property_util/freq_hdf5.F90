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

! Write the harmonic frequency analysis to the HDF5 file FREQH5.
!
! The first eight arguments are those of Freq_Molden, and the data they stand for are written in the same order, so that
! the two files can be compared row by row. The others are what only the caller knows:
!   Path       ANALYTIC, NUMERICAL or UPDATED (the optimizer's approximate Hessian)
!   Energy     energy at the geometry of the frequencies, written only if Do_Energy
!   nsRot      rotational symmetry number used by the thermochemistry, 0 if it did not run
!   Do_Intens  whether Intens holds computed intensities
! The file is created anew at every call.

subroutine Freq_HDF5(Freq,nFreq,Vectors,nVectors,nSym,Intens,mDisp,RedMas,Path,Energy,Do_Energy,nsRot,Do_Intens)

#ifdef _HDF5_
use mh5, only: mh5_close_dset, mh5_close_file, mh5_create_dset_int, mh5_create_dset_real, mh5_create_dset_str, mh5_create_file, &
               mh5_init_attr, mh5_put_dset
use Symmetry_Info, only: ChOper, iOper, nIrrep, Symmetry_Info_Get
use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nFreq, nVectors, nSym, mDisp(nSym), nsRot
real(kind=wp), intent(in) :: Freq(nFreq), Vectors(nVectors), Intens(nFreq), RedMas(nFreq), Energy
character(len=*), intent(in) :: Path
logical(kind=iwp), intent(in) :: Do_Energy, Do_Intens
#ifdef _HDF5_
integer(kind=iwp) :: dsetid, fileid, iAll, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iSpin, iStab(0:7), iUnique, lcDmp, MaxDCR, &
                     nAll_Atoms, nCnttp, nCoSet, nGen, nStab, nTmp, nUnique_Atoms
logical(kind=iwp) :: Found, Have_AtoB, Have_cDmp, Have_Iso
integer(kind=iwp), allocatable :: AtNum(:), AtNum_All(:), AtoB(:)
real(kind=wp), allocatable :: Charge(:), Coord(:,:), Coord_All(:,:), Mass(:), NMode(:,:,:)
character(len=8) :: Method
character(len=LenIn), allocatable :: Label(:)
character(len=LenIn+4), allocatable :: Label_All(:)
character(len=80), allocatable :: Bsl(:), Bsl_Old(:)
character(len=2*80), allocatable :: cDmp(:)
integer(kind=iwp), external :: iChxyz

call Symmetry_Info_Get()

call Get_iScalar('Unique atoms',nUnique_Atoms)
call Get_nAtoms_All(nAll_Atoms)

fileid = mh5_create_file('FREQH5')

if (Path == 'ANALYTIC') then
  call mh5_init_attr(fileid,'MOLCAS_MODULE','MCLR')
else
  call mh5_init_attr(fileid,'MOLCAS_MODULE','SLAPAF')
end if
call mh5_init_attr(fileid,'NSYM',nIrrep)
call mh5_init_attr(fileid,'NATOMS_UNIQUE',nUnique_Atoms)
if (nIrrep > 1) call mh5_init_attr(fileid,'NATOMS_ALL',nAll_Atoms)
call mh5_init_attr(fileid,'NFREQ',nFreq)
call mh5_init_attr(fileid,'HESSIAN_SOURCE',trim(Path))
! LDV: writing this file must never stop a calculation. The records below that the frequency analysis may not have needed
! before this point ('Multiplicity', 'Relax Method', 'Atom -> Basis', 'Isotopes', 'cDmp') are tested first, and what they feed
! is left out of the file when one is missing. The records that Freq_Molden reads at the same call need no test, and neither
! does 'Unique Atom Names', which MCLR and SLAPAF read when they start.
call Qpg_iScalar('Multiplicity',Found)
if (Found) then
  call Get_iScalar('Multiplicity',iSpin)
  call mh5_init_attr(fileid,'SPINMULT',iSpin)
end if
call Qpg_cArray('Relax Method',Found,nTmp)
if (Found) then
  call Get_cArray('Relax Method',Method,8)
  call mh5_init_attr(fileid,'RELAX_METHOD',trim(Method))
end if
! LDV: nsRot is 0 when the thermochemistry did not run; it is never 0 after Rotation has run.
if (nsRot > 0) call mh5_init_attr(fileid,'ROT_SYMMETRY_NUMBER',nsRot)
!                                                                      *
!***********************************************************************
!                                                                      *
! Symmetry-unique centers, as in the wavefunction files.
! LDV: pseudo-charge centers are not among 'Unique atoms', hence not among the centers written here, and neither are they
! among the atoms of freq.molden.

call mma_allocate(Label,nUnique_Atoms,label='Label')
call Get_cArray('Unique Atom Names',Label,LenIn*nUnique_Atoms)
dsetid = mh5_create_dset_str(fileid,'CENTER_LABELS',1,[nUnique_Atoms],LenIn)
call mh5_init_attr(dsetid,'DESCRIPTION','Unique center labels arranged as one [NATOMS_UNIQUE] block; '// &
                   'pseudo-charge centers are not included, as they are not atoms of the frequency analysis')
call mh5_put_dset(dsetid,Label)
call mh5_close_dset(dsetid)

! LDV: the RunFile holds the atomic numbers of the unique centers as reals, 0 for a dummy center.
call mma_allocate(Charge,nUnique_Atoms,label='Charge')
call Get_dArray('Nuclear charge',Charge,nUnique_Atoms)
call mma_allocate(AtNum,nUnique_Atoms,label='AtNum')
AtNum(:) = nint(Charge)
call mma_deallocate(Charge)
dsetid = mh5_create_dset_int(fileid,'CENTER_ATNUMS',1,[nUnique_Atoms])
call mh5_init_attr(dsetid,'DESCRIPTION','Atomic numbers of the unique centers, stored as array of size [NATOMS_UNIQUE]; '// &
                   '0 marks a dummy center')
call mh5_put_dset(dsetid,AtNum)
call mh5_close_dset(dsetid)

! LDV: the effective nuclear charge differs from the atomic number for a ghost atom of a real element (0), and for a center with
! an effective core potential; it is the only way to tell such a ghost atom from the file.
call Qpg_dArray('Effective nuclear Charge',Found,nTmp)
if (Found) then
  call mma_allocate(Charge,nUnique_Atoms,label='Charge')
  call Get_dArray('Effective nuclear Charge',Charge,nUnique_Atoms)
  dsetid = mh5_create_dset_real(fileid,'CENTER_CHARGES',1,[nUnique_Atoms])
  call mh5_init_attr(dsetid,'DESCRIPTION','Effective nuclear charges of the unique centers, stored as array of size '// &
                     '[NATOMS_UNIQUE]; it is 0 for a ghost atom of a real element, which keeps its atomic number in '// &
                     'CENTER_ATNUMS, and it is the charge seen by the electrons, less than the atomic number, for a center '// &
                     'with an effective core potential')
  call mh5_put_dset(dsetid,Charge)
  call mh5_close_dset(dsetid)
  call mma_deallocate(Charge)
end if

call mma_allocate(Coord,3,nUnique_Atoms,label='Coord')
call Get_dArray('Unique Coordinates',Coord,3*nUnique_Atoms)
dsetid = mh5_create_dset_real(fileid,'CENTER_COORDINATES',2,[3,nUnique_Atoms])
call mh5_init_attr(dsetid,'DESCRIPTION','Unique center coordinates in bohr, matrix of size [NATOMS_UNIQUE,3], stored with '// &
                   'atom index varying slowest; the geometry the frequencies belong to')
call mh5_put_dset(dsetid,Coord)
call mh5_close_dset(dsetid)

! LDV: the masses are those of the RunFile, which both frequency analyses use for the real atoms; a dummy center has mass 0.
call Qpg_iArray('Atom -> Basis',Have_AtoB,nTmp)
call Qpg_dArray('Isotopes',Have_Iso,nTmp)
if (Have_AtoB .and. Have_Iso) then
  call mma_allocate(Mass,nUnique_Atoms,label='Mass')
  call Get_Mass(Mass,nUnique_Atoms)
  dsetid = mh5_create_dset_real(fileid,'CENTER_MASSES',1,[nUnique_Atoms])
  call mh5_init_attr(dsetid,'DESCRIPTION','Nuclear masses in electron masses (atomic units of mass), stored as array of '// &
                     'size [NATOMS_UNIQUE], including any isotope set in the input; 0 marks a dummy center; the images of '// &
                     'a unique center share its mass')
  call mh5_put_dset(dsetid,Mass)
  call mh5_close_dset(dsetid)
  call mma_deallocate(Mass)
end if

! LDV: 'Atom -> Basis' gives the centre type of each unique atom, auxiliary, pseudo-charge and fragment types being left out
! of it, and 'cDmp' holds the two labels of each centre type, 80 characters each, the library one then the one in the input.
call Qpg_cArray('cDmp',Have_cDmp,lcDmp)
if (Have_AtoB .and. Have_cDmp) then
  call mma_allocate(AtoB,nUnique_Atoms,label='AtoB')
  call Get_iArray('Atom -> Basis',AtoB,nUnique_Atoms)
  nCnttp = lcDmp/(2*80)
  call mma_allocate(cDmp,nCnttp,label='cDmp')
  call Get_cArray('cDmp',cDmp,lcDmp)
  call mma_allocate(Bsl,nUnique_Atoms,label='Bsl')
  call mma_allocate(Bsl_Old,nUnique_Atoms,label='Bsl_Old')
  do iUnique=1,nUnique_Atoms
    Bsl(iUnique) = cDmp(AtoB(iUnique))(1:80)
    Bsl_Old(iUnique) = cDmp(AtoB(iUnique))(81:160)
  end do
  call mma_deallocate(cDmp)
  call mma_deallocate(AtoB)

  dsetid = mh5_create_dset_str(fileid,'CENTER_BASIS_LABELS',1,[nUnique_Atoms],80)
  call mh5_init_attr(dsetid,'DESCRIPTION','Basis set labels of the unique centers, arranged as one [NATOMS_UNIQUE] block of '// &
                     '80 characters, stored as kept by the program: the library form that the SEWARD output prints as the '// &
                     'basis set label, in upper case, with the library name as translated by the basis set table, ending in '// &
                     'a run of dots and then blanks; the images of a unique center share its basis set')
  call mh5_put_dset(dsetid,Bsl)
  call mh5_close_dset(dsetid)

  dsetid = mh5_create_dset_str(fileid,'CENTER_BASIS_INPUT_LABELS',1,[nUnique_Atoms],80)
  call mh5_init_attr(dsetid,'DESCRIPTION','Basis set labels of the unique centers as given in the input, arranged as one '// &
                     '[NATOMS_UNIQUE] block of 80 characters, stored as kept by the program (the label is padded with dots '// &
                     'to 80 characters), in the form that the user wrote, which may differ from CENTER_BASIS_LABELS in the '// &
                     'name of the basis set')
  call mh5_put_dset(dsetid,Bsl_Old)
  call mh5_close_dset(dsetid)
  call mma_deallocate(Bsl)
  call mma_deallocate(Bsl_Old)
end if
!                                                                   *
!***********************************************************************
!                                                                      *
! All centers, only with symmetry, in the order of the coordinates and normal modes of Freq_Molden: each unique center is
! followed by its images, in the order of its coset (the loop is that of Get_Coord_All_).

if (nIrrep > 1) then
  nGen = 0
  if (nIrrep == 2) nGen = 1
  if (nIrrep == 4) nGen = 2
  if (nIrrep == 8) nGen = 3
  if (nGen >= 1) iGen(1) = iOper(1)
  if (nGen >= 2) iGen(2) = iOper(2)
  if (nGen == 3) iGen(3) = iOper(4)

  call mma_allocate(Label_All,nAll_Atoms,label='Label_All')
  call mma_allocate(AtNum_All,nAll_Atoms,label='AtNum_All')
  iAll = 0
  MaxDCR = 0
  do iUnique=1,nUnique_Atoms
    iChAtom = iChxyz(Coord(:,iUnique),iGen,nGen)
    call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
    nCoSet = nIrrep/nStab
    do iCo=0,nCoSet-1
      iAll = iAll+1
      Label_All(iAll) = Label(iUnique)//':'//ChOper(iCoSet(iCo,0))
      AtNum_All(iAll) = AtNum(iUnique)
    end do
  end do

  dsetid = mh5_create_dset_str(fileid,'DESYM_CENTER_LABELS',1,[nAll_Atoms],LenIn+4)
  call mh5_init_attr(dsetid,'DESCRIPTION','Desymmetrized center labels arranged as one [NATOMS_ALL] block, each unique '// &
                     "center followed by its images, labelled with the symmetry operation ('E', 'x', 'xy', ...) that "// &
                     'generates them')
  call mh5_put_dset(dsetid,Label_All)
  call mh5_close_dset(dsetid)

  dsetid = mh5_create_dset_int(fileid,'DESYM_CENTER_ATNUMS',1,[nAll_Atoms])
  call mh5_init_attr(dsetid,'DESCRIPTION','Desymmetrized atomic numbers, stored as array of size [NATOMS_ALL]; 0 marks a '// &
                     'dummy center')
  call mh5_put_dset(dsetid,AtNum_All)
  call mh5_close_dset(dsetid)

  call mma_deallocate(Label_All)
  call mma_deallocate(AtNum_All)
end if
call mma_deallocate(Label)
call mma_deallocate(AtNum)
call mma_deallocate(Coord)

! All coordinates, as written to freq.molden (they are the unique ones when there is no symmetry).
if (nIrrep > 1) then
  call mma_allocate(Coord_All,3,nAll_Atoms,label='Coord_All')
  call Get_Coord_All(Coord_All,nAll_Atoms)
  dsetid = mh5_create_dset_real(fileid,'DESYM_CENTER_COORDINATES',2,[3,nAll_Atoms])
  call mh5_init_attr(dsetid,'DESCRIPTION','Desymmetrized coordinates in bohr, size [NATOMS_ALL,3], stored with atom index '// &
                     'varying slowest; the geometry the frequencies belong to')
  call mh5_put_dset(dsetid,Coord_All)
  call mh5_close_dset(dsetid)
  call mma_deallocate(Coord_All)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Frequencies, normal modes and reduced masses, as in freq.molden.

dsetid = mh5_create_dset_real(fileid,'FREQUENCIES',1,[nFreq])
call mh5_init_attr(dsetid,'DESCRIPTION','Harmonic frequencies in cm-1, vector of size [NFREQ] in the order of '// &
                   'freq.molden; a negative value is an imaginary frequency; the root attributes are NFREQ, the number '// &
                   'of frequencies; NSYM, the number of operations of the point group used in the calculation; '// &
                   'NATOMS_UNIQUE and, only when NSYM is larger than 1, NATOMS_ALL, the numbers of symmetry-unique and '// &
                   'of all atoms; MOLCAS_MODULE, the module that wrote the file, MCLR for an analytic Hessian and SLAPAF '// &
                   'otherwise; MOLCAS_VERSION, the version of the code; and HESSIAN_SOURCE, ANALYTIC for an analytic '// &
                   'Hessian, NUMERICAL for a Hessian obtained by finite differences of gradients, and UPDATED for the '// &
                   'approximate Hessian updated by the geometry optimizer, which is not a computed second derivative '// &
                   'and must not be used for thermochemistry; SPINMULT, the spin multiplicity, and RELAX_METHOD, the '// &
                   'label of the method that produced the wavefunction, each absent if the calculation did not record it; '// &
                   'and ROT_SYMMETRY_NUMBER, present only when thermochemistry was done, the rotational symmetry number '// &
                   'that it used, which is the value given in the input, or 1 if none was given, except that '// &
                   'thermochemistry uses 2 for a homonuclear diatomic molecule whatever the input says, and is not '// &
                   'determined from the geometry')
call mh5_put_dset(dsetid,Freq)
call mh5_close_dset(dsetid)

call mma_allocate(NMode,3,nAll_Atoms,nFreq,label='NMode')
NMode(:,:,:) = Zero
call Get_NMode_All(Vectors,nVectors,nFreq,nUnique_Atoms,NMode,nAll_Atoms,mDisp)
dsetid = mh5_create_dset_real(fileid,'NORMAL_MODES',3,[3,nAll_Atoms,nFreq])
call mh5_init_attr(dsetid,'DESCRIPTION','Normal modes as written to freq.molden, array of size [NFREQ,NATOMS,3], stored '// &
                   'with mode index varying slowest, then atom index, where NATOMS is NATOMS_ALL when NSYM is larger '// &
                   'than 1 and NATOMS_UNIQUE otherwise, and the atoms are in the order of DESYM_CENTER_COORDINATES, or '// &
                   'of CENTER_COORDINATES without symmetry; the normalization depends on the path: for an analytic Hessian '// &
                   'each mode has the sum over all atoms of |v|^2 equal to 1, and for a numerical Hessian or the '// &
                   "optimizer's updated Hessian with the default coordinates of the optimizer it has the sum over all atoms "// &
                   'of m|v|^2 equal to 1 with m the mass in u, while other coordinate options can give another '// &
                   'normalization, so a tool that depends on a convention should normalize the modes itself')
call mh5_put_dset(dsetid,NMode)
call mh5_close_dset(dsetid)
call mma_deallocate(NMode)

dsetid = mh5_create_dset_real(fileid,'REDUCED_MASSES',1,[nFreq])
call mh5_init_attr(dsetid,'DESCRIPTION','Reduced masses of the normal modes in u, vector of size [NFREQ] in the order of '// &
                   'FREQUENCIES; each is the sum over all atoms of m|v|^2 divided by the sum over all atoms of |v|^2, where '// &
                   'v is the displacement of the atom in the mode as stored in NORMAL_MODES and m the mass of the atom in u '// &
                   '(the mass in CENTER_MASSES converted from electron masses), so that it does not depend on how the '// &
                   'mode is normalized')
call mh5_put_dset(dsetid,RedMas)
call mh5_close_dset(dsetid)

! LDV: only written when the intensities were computed, never as zeros standing for 'not computed'.
if (Do_Intens) then
  dsetid = mh5_create_dset_real(fileid,'IR_INTENSITIES',1,[nFreq])
  call mh5_init_attr(dsetid,'DESCRIPTION','Infrared intensities of the normal modes in km/mol, vector of size [NFREQ] '// &
                     'in the order of FREQUENCIES; zero for a mode of an irreducible representation without a dipole '// &
                     'component; present only when the intensities were computed')
  call mh5_put_dset(dsetid,Intens)
  call mh5_close_dset(dsetid)
end if

! LDV: the energy at the geometry of the frequencies, for the relaxation root; the caller says when it is not available.
if (Do_Energy) then
  dsetid = mh5_create_dset_real(fileid,'ENERGY')
  call mh5_init_attr(dsetid,'DESCRIPTION','Total energy in hartree at the geometry of the frequencies, for the '// &
                     'relaxation root, the root whose gradient and Hessian are computed; absent when the calculation '// &
                     'averages two states')
  call mh5_put_dset(dsetid,Energy)
  call mh5_close_dset(dsetid)
end if

call mh5_close_file(fileid)
#else
#include "macros.fh"
unused_var(Freq)
unused_var(nFreq)
unused_var(Vectors)
unused_var(nVectors)
unused_var(nSym)
unused_var(Intens)
unused_var(mDisp)
unused_var(RedMas)
unused_var(Path)
unused_var(Energy)
unused_var(Do_Energy)
unused_var(nsRot)
unused_var(Do_Intens)
#endif

end subroutine Freq_HDF5
