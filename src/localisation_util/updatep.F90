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
! Copyright (C) Yannick Carissan                                       *
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine UpdateP(PACol,Name,nBas_Start,nOrb2Loc,nAtoms,PA,gamma_rot,iMO_s,iMO_t,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Reduce operation count and use BLAS.

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
integer nBas_Start(*)
real*8 PACol(nOrb2Loc,2)
real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
character*(LENIN8) Name(*), PALbl
logical Debug

cosg = cos(gamma_rot)
sing = sin(gamma_rot)
cos2g = cosg*cosg
sin2g = sing*sing
cosing = cosg*sing

do iAt=1,nAtoms
  !call RecPrt('PA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)

  ! Copy out the PAss, PAtt, and PAst elements.

  pass = PA(iMO_s,iMO_s,iAt)
  PAst = PA(iMO_s,iMO_t,iAt)
  PAtt = PA(iMO_t,iMO_t,iAt)
  !write(6,*) 'updateP:',PAss,PAst,PAtt
# if defined (_DEBUGPRINT_)
  PAts = PA(iMO_t,iMO_s,iAt)
  Tst = PAst-PAts
  if (abs(Tst) > 1.0d-14) then
    write(6,*) 'Broken symmetry in UpdateP!!'
    write(6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(6,*) 'PAst = ',PAst
    write(6,*) 'PAts = ',PAts
    write(6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateP','Broken symmetry!',' ')
  end if
# endif

  ! Copy out columns s and t of PA.

  call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PACol(1,1),1)
  call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PACol(1,2),1)

  ! Compute transformed columns.

  call dScal_(nOrb2Loc,cosg,PA(1,iMO_s,iAt),1)
  call dAXPY_(nOrb2Loc,sing,PACol(1,2),1,PA(1,iMO_s,iAt),1)
  call dScal_(nOrb2Loc,cosg,PA(1,iMO_t,iAt),1)
  call dAXPY_(nOrb2Loc,-sing,PACol(1,1),1,PA(1,iMO_t,iAt),1)

  ! Compute PAss, PAtt, PAst, and PAts (= PAst).

  PA(iMO_s,iMO_s,iAt) = pass*cos2g+PAtt*sin2g+Two*PAst*cosing
  PA(iMO_t,iMO_s,iAt) = (PAtt-pass)*cosing+PAst*(cos2g-sin2g)
  PA(iMO_s,iMO_t,iAt) = PA(iMO_t,iMO_s,iAt)
  PA(iMO_t,iMO_t,iAt) = PAtt*cos2g+pass*sin2g-Two*PAst*cosing

  ! Copy columns to rows.

  call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PA(iMO_s,1,iAt),nOrb2Loc)
  call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PA(iMO_t,1,iAt),nOrb2Loc)

end do

if (Debug) then
  write(6,*) 'In UpdateP'
  write(6,*) '----------'
  do iAt=1,nAtoms
    PALbl = 'PA__'//Name(nBas_Start(iAt))(1:LENIN)
    call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
  end do
end if

return

end subroutine UpdateP
