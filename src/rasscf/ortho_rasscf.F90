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
! Copyright (C) 1992, Per Ake Malmqvist                                *
!               1992, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Ortho_RASSCF(SMAT,SCRATCH,CMO,TEMP)
!***********************************************************************
!                                                                      *
!     purpose: Orthogonalize MOs (one symmetry block at a time)        *
!                                                                      *
!     calling arguments:                                               *
!     Smat    : overlap matrix                                         *
!     CMO     : MO-coefficients                                        *
!     Temp    : temporary work space                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.-AA. Malmqvist and M.P. Fuelscher                              *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use general_data, only: LOWDIN_ON, NBAS, NDEL, NSYM
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: Smat(*), SCRATCH(*), CMO(*), Temp(*)
integer(kind=iwp) :: i_Component, i_Opt, i_RC, i_SymLbl, iBas, iOcc, ip_CMO, ip_SMat, iSym
character(len=8) :: Label
#include "warnings.h"

!                                                                      *
!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
! Read overlap matrix SMAT:

i_Rc = 0
i_Opt = ibset(ibset(0,sNoOri),sNoNuc)
i_Component = 1
i_SymLbl = 1
Label = 'Mltpl  0'
call RdOne(i_Rc,i_Opt,Label,i_Component,Smat,i_SymLbl)
if (i_Rc /= 0) then
  write(u6,*) ' ORTHO could not read overlaps from ONEINT.'
  write(u6,*) ' RASSCF is trying to orthonormalize orbitals but'
  write(u6,*) ' could not read overlaps from ONEINT. Something'
  write(u6,*) ' is wrong with the file, or possibly with the'
  write(u6,*) ' program. Please check.'
  call Quit(_RC_IO_ERROR_READ_)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Orthonormalize symmetry blocks:

ip_Smat = 1
ip_CMO = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  iOcc = iBas-nDel(iSym)
  if (iBas > 0) then

    call SQUARE(SMAT(ip_Smat),Temp,1,iBas,iBas)

    !call RecPrt('S',' ',Temp,iBas,iBas)
    !call RecPrt('CMO',' ',CMO(ip_CMO),iBas,iBas)

    if (Lowdin_ON) then

      ! compute C^T*S*C = W  (overlap in MO basis)

      call DGEMM_('T','N',iOcc,iBas,iBas,One,CMO(ip_CMO),iBas,Temp,iBas,Zero,SCRATCH,iOcc)
      call DGEMM_('N','N',iOcc,iOcc,iBas,One,SCRATCH,iOcc,CMO(ip_CMO),iBas,Zero,Temp,iOcc)

      ! compute W^-1/2

      call Lowdin_LP(Temp,SCRATCH,iOcc)

      ! compute C' = C*W^-1/2

      call DGEMM_('N','N',iBas,iOcc,iOcc,One,CMO(ip_CMO),iBas,SCRATCH,iOcc,Zero,Temp,iBas)

      ! PAM March 2016: Probable bugfix needed (Thanks, Liviu!)
      ! not affecting any tests (!)
      ! by adding the following line:
      call DCOPY_(iBas*iOcc,Temp,1,CMO(ip_CMO),1)
    else

      call ORTHO1(Temp,CMO(ip_CMO),SCRATCH,iBas,iOcc)
    end if

    !call RecPrt('CMO',' ',CMO(ip_CMO),iBas,iBas)

    ip_Smat = ip_Smat+nTri_Elem(iBas)
    ip_CMO = ip_CMO+iBas*iBas
  end if
end do

end subroutine Ortho_RASSCF
