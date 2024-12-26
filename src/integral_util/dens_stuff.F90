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
! This should probably be changes to arrays and pointers at some later point.
Module Dens_stuff
use Definitions, only: iwp
Private
integer(kind=iwp), public, target :: mDCRij=1,mDCRkl=1,mDCRik=1,mDCRil=1,mDCRjk=1,mDCRjl=1
integer(kind=iwp), public, target ::  ipDij, ipDkl, ipDik, ipDil, ipDjk, ipDjl
integer(kind=iwp), public, target :: ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl
integer(kind=iwp), public, target ::   mDij,  mDkl,  mDik,  mDil,  mDjk,  mDjl
End Module Dens_stuff
