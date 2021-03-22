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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,nFro,nOcc,nVir,nDel,CMOinp,CMOout)

implicit none
integer nSym
integer nBas(nSym)
integer nOrb(nSym)
integer nFro(nSym)
integer nOcc(nSym)
integer nVir(nSym)
integer nDel(nSym)
real*8 CMOinp(*)
real*8 CMOout(*)

integer iSym, ip1, ip2, l

ip1 = 1
ip2 = 1
do iSym=1,nSym
  l = nBas(iSym)*nOrb(iSym)
  call dCopy_(l,CMOinp(ip1),1,CMOout(ip2),1)
  call fZero(CMOout(ip2+l),nBas(iSym)*nBas(iSym)-l)
  ip1 = ip1+l
  ip2 = ip2+nBas(iSym)*nBas(iSym)
end do

return

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(nDel)
  call Unused_integer_array(nFro)
  call Unused_integer_array(nOcc)
  call Unused_integer_array(nVir)
end if

end subroutine ChoRPA_MOTra_ReorderCMO
