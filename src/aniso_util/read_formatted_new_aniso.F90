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

subroutine read_formatted_new_aniso(input_file_name,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO,eso_au, &
                                    esfs_au)

use Constants, only: auTocm, gElectron
use Definitions, only: u6

implicit none
integer, intent(inout) :: nss, nstate
integer, intent(out) :: multiplicity(nstate)
real(kind=8), intent(out) :: eso(nss), esfs(nstate)
real(kind=8), intent(out) :: eso_au(nss), esfs_au(nstate)
real(kind=8), intent(out) :: edmom(3,nstate,nstate)
real(kind=8), intent(out) :: angmom(3,nstate,nstate)
real(kind=8), intent(out) :: amfi(3,nstate,nstate)
complex(kind=8), intent(out) :: MM(3,nss,nss)
complex(kind=8), intent(out) :: MS(3,nss,nss)
complex(kind=8), intent(out) :: ML(3,nss,nss)
! electric dipole moment
complex(kind=8), intent(out) :: DM(3,nss,nss)
complex(kind=8), intent(out) :: U(nss,nss)
complex(kind=8), intent(out) :: HSO(nss,nss)
character(len=180) :: input_file_name
! local variables:
integer :: LuAniso, IsFreeUnit
real(kind=8), parameter :: g_e = -gElectron
external :: IsFreeUnit
logical :: DBG

dbg = .false.
if (dbg) write(u6,'(A)') 'Entering read_formatted_aniso_new'
! read the data file:
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,input_file_name)
call read_magnetic_moment(LuAniso,nss,MM,dbg)
call read_electric_moment(LuAniso,nss,DM,dbg)
call read_spin_moment(LuAniso,nss,MS,dbg)
call read_angmom(LuAniso,nstate,angmom,dbg)
call read_amfi(LuAniso,nstate,amfi,dbg)
call read_edipmom(LuAniso,nstate,EDMOM,dbg)
call read_multiplicity(LuAniso,nstate,multiplicity,dbg)
call read_eso(LuAniso,nss,eso_au,dbg)
call read_esfs(LuAniso,nstate,esfs_au,dbg)
call read_hso(LuAniso,nss,HSO,dbg)
call read_eigen(LuAniso,nss,U,dbg)

! compute the relative spin-orbit energies in cm-1
eso(:) = (eso_au(:)-eso_au(1))*auTocm

! compute the relative spin-free energies in cm-1
esfs(:) = (esfs_au(:)-esfs_au(1))*auTocm

! compute the orbital moment
ML(:,:,:) = -MM(:,:,:)-g_e*MS(:,:,:)

!read_nmult,
!read_imult,
!read_format,
!read_nroot,
!read_szproj,

close(LuAniso)

return

end subroutine read_formatted_new_aniso
