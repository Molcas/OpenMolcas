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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               1995,1996, Pavel Neogrady                              *
!***********************************************************************

subroutine mod2(nsym,nish,nash,nssh,norb,fi,eps)
! this routine defines fi(p,q) = delta(p,q).eps(p)
! and redefines nish=nish+nash, nash=0
!
! this is suitable for closed shell case

integer nsym
integer nish(1:8)
integer nash(1:8)
integer nssh(1:8)
integer norb(1:8)
real*8 fi(*)
real*8 eps(*)
! help variables
integer p, q, isym, pq, padd

!1 redefine foki

pq = 0
padd = 0
do isym=1,nsym

  do p=1,norb(isym)
    do q=1,p
      pq = pq+1
      if (p == q) then
        fi(pq) = eps(padd+p)
      else
        fi(pq) = 0.0d0
      end if
    end do
  end do

  padd = padd+norb(isym)
end do

!2 redefine n's

do isym=1,nsym
  nish(isym) = nish(isym)+nash(isym)
  nash(isym) = 0
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(nssh)

end subroutine mod2
