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

subroutine write_new_formatted_aniso(nss,nstate,multiplicity,eso_au,esfs_au,U,MM,MS,DM,angmom,edmom,amfi,HSO)

use Constants, only: Angstrom

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: nss, nstate, multiplicity(nstate)
real(kind=8), intent(in) :: eso_au(nss), esfs_au(nstate)
real(kind=8), intent(in) :: angmom(3,nstate,nstate)
real(kind=8), intent(in) :: edmom(3,nstate,nstate)
real(kind=8), intent(in) :: amfi(3,nstate,nstate)
complex(kind=8), intent(in) :: MM(3,nss,nss)
complex(kind=8), intent(in) :: MS(3,nss,nss)
complex(kind=8), intent(in) :: DM(3,nss,nss)
complex(kind=8), intent(in) :: U(nss,nss)
complex(kind=8), intent(in) :: HSO(nss,nss)
! local stuff
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
integer :: njob, mxjob, mult, iss, ipar, ist
integer :: data_file_format
integer :: i, IsFreeUnit, Lu, Lutmp
!character(len=30) :: fmt_int, fmt_real, fmt_key
external :: IsFreeUnit
integer, allocatable :: szproj(:), jbnum(:), mltplt(:)
integer, allocatable :: nroot(:)
character(len=128) :: Filename
character(len=1024) :: molcas, fname, molcasversion
logical :: dbg
integer :: nAtoms, iAt, l
character(len=LENIN) :: AtomLbl(MxAtom)
real*8, allocatable :: xyz(:,:)

dbg = .false.

!-----------------------------------------------------------------------
! some preparations
njob = 0
mxjob = 0
call get_iScalar('NJOB_SINGLE',njob)
call get_iScalar('MXJOB_SINGLE',mxjob)
! allocate temporary memory:
call mma_allocate(jbnum,nstate,'jbnum')
call mma_allocate(mltplt,mxjob,'mltplt')
! get the information from RUNFILE:
mltplt = 0
jbnum = 0
call get_iArray('MLTP_SINGLE',mltplt,mxjob)
call get_iArray('JBNUM_SINGLE',jbnum,nstate)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! prepare the szproj index table
call mma_allocate(szproj,nss,'szproj')
szproj(1:nss) = 0
iss = 0
ipar = mod(multiplicity(1),2)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) go to 310
    Iss = Iss+1
    szproj(iss) = I
310 continue
  end do ! i
end do ! ist

call get_iScalar('MXJOB_SINGLE',mxjob)
call mma_allocate(nroot,mxjob,'nroot')
nroot = 0
call get_iArray('NSTAT_SINGLE',nroot,mxjob)

!-----------------------------------------------------------------------
! Get the MOLCAS version: index table
call getenvf('MOLCAS ',molcas)
write(fname,'(A)') trim(molcas)//'/.molcasversion'
Lutmp = IsFreeUnit(89)
call molcas_open(Lutmp,fname)
read(Lutmp,'(A180)') molcasversion
close(Lutmp)

!-----------------------------------------------------------------------
! add coordinates
nAtoms = 0
call Get_iScalar('Unique atoms',nAtoms)
call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
call mma_allocate(xyz,3,8*nAtoms)
call Get_dArray('Unique Coordinates',xyz,3*nAtoms)

!-----------------------------------------------------------------------
! write the data to the new aniso file
FileName = 'ANISOFILE'
Lu = IsFreeUnit(81)

call molcas_open(Lu,FileName)
data_file_format = 2021

!fmt_key = '(A)'
!fmt_real = '(5ES22.14,1x)'
!fmt_int = '(40(I0,1x))'

write(Lu,'(        A)') '# OPENMOLCAS interface to ANISO'
!-----------------------------------------------------------------------
! ORIGIN of DATA file
write(Lu,'(        A)') '$source '
write(Lu,'(       2A)') 'MOLCAS  ',trim(molcas)
write(Lu,'(       2A)') 'VERSION ',trim(molcasversion)
write(Lu,'(        A)')
!-----------------------------------------------------------------------
! DATA FILE FORMAT VERSION:
write(Lu,'(        A)') '$format '
write(Lu,'(40(I0,1x))') data_file_format
write(Lu,'(        A)')
!-----------------------------------------------------------------------
! Atom labels and coordinates
write(Lu,'(        A)') '$natoms '
write(Lu,'(40(I0,1x))') nAtoms
write(Lu,'(        A)')
write(Lu,'(        A)') '$atomlbl'
write(Lu,'(40(I0,1x))') nAtoms
write(Lu,'(40(A8,1x))') (AtomLbl(iAt),iAt=1,nAtoms)
write(Lu,'(        A)')
write(Lu,'(        A)') '$coords (in Angstrom)'
write(Lu,'(40(I0,1x))') nAtoms
do iAt=1,nAtoms
  write(Lu,'(i3,1x,A8,1x,3(ES22.14,1x))') iAt,AtomLbl(iAt),(Angstrom*XYZ(l,iAt),l=1,3)
end do
write(Lu,'(        A)')
!-----------------------------------------------------------------------
! Number of spin orbit states
call write_nss(LU,nss,dbg)
!-----------------------------------------------------------------------
! Number of spin free states
call write_nstate(LU,nstate,dbg)
!-----------------------------------------------------------------------
! Number of spin multiplicities
call write_nmult(LU,njob,dbg)
!-----------------------------------------------------------------------
! Values of the spin multiplicity
call write_imult(LU,njob,mltplt,dbg)
!-----------------------------------------------------------------------
! Number of roots in each spin multiplicity
call write_nroot(LU,njob,nroot,dbg)
!-----------------------------------------------------------------------
! SZ projection of the states defining the ording
call write_szproj(LU,nss,szproj,dbg)
!-----------------------------------------------------------------------
! spin multiplicity for all states
call write_multiplicity(LU,nstate,multiplicity,dbg)
!-----------------------------------------------------------------------
! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
! atomic units
call write_eso(LU,nss,eso_au,dbg)
!-----------------------------------------------------------------------
! Eigenvalues of the HBO matrix (spin-free energies from RASSCF)
! atomic units
call write_esfs(LU,nstate,esfs_au,dbg)
!-----------------------------------------------------------------------
! Angular momentum operator in the basis of spin free states,
! component X,Y,Z
call write_angmom(LU,nstate,angmom,dbg)
!-----------------------------------------------------------------------
! AMFI/SOMF operator in the basis of spin free states
! (see eq. 35 in:
! D. Ganyushin and F. Neese, J. Chem. Phys. 138, 104113, 2013).
! Note that the meaning of AMFI integrals is different in
! MOLCAS and ORCA
call write_amfi(LU,nstate,amfi,dbg)
!-----------------------------------------------------------------------
! Electric transition-dipole moments in the basis of spin free states
! The dipole moments, where bra and ket have the same state, are not computed. ??
call write_edipmom(LU,nstate,edmom,dbg)
!-----------------------------------------------------------------------
! Magnetic dipole moment in the basis of SO states
call write_magnetic_moment(LU,nss,MM,dbg)
!-----------------------------------------------------------------------
! Spin moment in the basis of SO states
call write_spin_moment(Lu,nss,MS,dbg)
!-----------------------------------------------------------------------
! Electric transition-dipole moments in the basis of SO states
call write_electric_moment(LU,nss,DM,dbg)
!-----------------------------------------------------------------------
! Eigenvectors of the HBO+SOC matrix
call write_eigen(LU,nss,U,dbg)
!-----------------------------------------------------------------------
! The HBO+SOC matrix, atomic units
call write_hso(LU,nss,HSO,dbg)
!-----------------------------------------------------------------------
flush(Lu)
close(Lu)

call mma_deallocate(szproj)
call mma_deallocate(jbnum)
call mma_deallocate(mltplt)
call mma_deallocate(nroot)
call mma_deallocate(xyz)

return

end subroutine write_new_formatted_aniso
