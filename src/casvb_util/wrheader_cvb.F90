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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: recn
integer(kind=iwp) :: norb1, nbas_mo1, nvb1, kbasiscvb1, ioffs_orbs, ioffs_cvb, ioffs_orbsao, ioffs_orbslao
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iheader
integer(kind=iwp), external :: ihlf_cvb, mstackiz_cvb
logical(kind=iwp), parameter :: debug = .false.

ioffs_orbs = ihlf_cvb(10)
ioffs_cvb = ioffs_orbs+norb1*norb1
ioffs_orbsao = ioffs_cvb+nvb1
ioffs_orbslao = ioffs_orbsao+norb1*nbas_mo1
!call reserv_cvb(ioffs_orbslao+norb1*nbas_mo1,recn)
iheader = mstackiz_cvb(10)
iwork(iheader) = norb1
iwork(1+iheader) = nbas_mo1
iwork(2+iheader) = nvb1
iwork(3+iheader) = kbasiscvb1
iwork(5+iheader) = ioffs_orbs
iwork(6+iheader) = ioffs_cvb
iwork(7+iheader) = ioffs_orbsao
iwork(8+iheader) = ioffs_orbslao
call wri_cvb(iwork(iheader),10,recn,0)
call mfreei_cvb(iheader)
if (debug) then
  write(u6,*) ' wrheader :'
  write(u6,*) ' ----------'
  write(u6,*) ' norb1         :',norb1
  write(u6,*) ' nbas_mo1      :',nbas_mo1
  write(u6,*) ' nvb1          :',nvb1
  write(u6,*) ' kbasiscvb1    :',kbasiscvb1
  write(u6,*) ' ioffs_orbs    :',ioffs_orbs
  write(u6,*) ' ioffs_cvb     :',ioffs_cvb
  write(u6,*) ' ioffs_orbsao  :',ioffs_orbsao
  write(u6,*) ' ioffs_orbslao :',ioffs_orbslao
end if

return

end subroutine wrheader_cvb
