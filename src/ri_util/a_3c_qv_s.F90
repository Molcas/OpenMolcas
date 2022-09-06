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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine A_3C_Qv_s(A_3C,Qv,Rv,nMuNu,nI,nK,QMode)
!***********************************************************************
!
!     Author:  F. Aquilante
!
!***********************************************************************

implicit real*8(a-h,o-z)
real*8 A_3C(nMuNu,*), Qv(nI,nK), Rv(nMuNu,*)
character*1 QMode

if (QMode == 'N') then

  call DGEMM_('N','N',nMuNu,nK,nI,1.0d0,A_3C,nMuNu,Qv,nI,0.0d0,Rv,nMuNu)

else if (QMode == 'T') then

  call DGEMM_('N','T',nMuNu,nI,nK,1.0d0,A_3C,nMuNu,Qv,nI,1.0d0,Rv,nMuNu)  ! note that Rv is accumulated

else

  call WarningMessage(2,'A_3C_Qv_s: illegal QMode!')
  call Abend()
end if

return

end subroutine A_3C_Qv_s
