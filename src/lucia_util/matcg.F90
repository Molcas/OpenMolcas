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

subroutine MATCG(CIN,COUT,NROWI,NROWO,IROWI1,NGCOL,IGAT,GATSGN)
! Gather columns of CIN with phase
!
! COUT(IR,IC) = GATSGN(IC)*CIN(IR+IROWI1-1,IGAT(IC)) if IGAT(IC) /= 0
! COUT(IR,IC) = 0                                    if IGAT(IC) /= 0

implicit real*8(A-H,O-Z)
integer IGAT(*)
dimension GATSGN(*)
dimension CIN(NROWI,*), COUT(NROWO,*)

!write(6,*) ' MATCG NROWI,NROWO,IROWI1,NGCOL'
!write(6,*) NROWI,NROWO,IROWI1,NGCOL
do IG=1,NGCOL
  !write(6,*) ' igat,sign ',IGAT(IG),GATSGN(IG)
  if (IGAT(IG) == 0) then
    do IR=1,NROWO
      COUT(IR,IG) = 0.0d0
    end do
  else
    IGFRM = IGAT(IG)
    SIGN = GATSGN(IG)
    do IR=1,NROWO
      COUT(IR,IG) = SIGN*CIN(IROWI1-1+IR,IGFRM)
    end do
  end if
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' Column gathered matrix'
  call WRTMAT(COUT,NROWO,NGCOL,NROWO,NGCOL)
end if

end subroutine MATCG
