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

integer function nBas_Eff(NrExp,NrBas,Exp,Cff,nExp_Eff)

implicit real*8(a-h,o-z)
real*8 exp(NrExp), Cff(NrExp,NrBas)

nBas_Eff = NrBas

do iBas=1,NrBas

  do iExp=1,nExp_Eff

    if (Cff(iExp,iBas) /= 0.0d0) then
      nBas_Eff = NrBas-iBas+1
      return
    end if

  end do

end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Exp)

end function nBas_Eff
