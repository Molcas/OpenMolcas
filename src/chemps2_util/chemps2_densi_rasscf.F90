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
! Copyright (C) 2016, Sebastian Wouters                                *
!               2016, Quan Phung                                       *
!***********************************************************************
! CheMPS2-Molcas main interface
! Based on Block interface, writen by N. Nakatani
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016

subroutine CHEMPS2_DENSI_RASSCF(jRoot,D,DS,PS,PA,PT)
implicit real*8(A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"

dimension D(NACPAR), DS(NACPAR), PS(NACPR2), PA(NACPR2)
dimension PT(NAC,NAC,NAC,NAC)

call DCOPY_(NACPAR,0.0d0,0,D,1)
call DCOPY_(NACPAR,0.0d0,0,DS,1)
call DCOPY_(NACPR2,0.0d0,0,PS,1)
call DCOPY_(NACPR2,0.0d0,0,PA,1)

if (NACTEL > 1) then
  call chemps2_load2pdm(NAC,PT,jRoot)
  IJ_pack = 1
  do J=1,NAC
    do I=1,J
      D1sum = 0.0d0
      do K=1,NAC
        D1sum = D1sum+PT(K,K,I,J)
      end do
      D(IJ_pack) = D1sum/(NACTEL-1)
      IJ_pack = IJ_pack+1
    end do
  end do

  IJKL_pack = 0
  do I=1,NAC
    do J=1,I
      do K=1,I
        LLIM = K
        if (K == I) LLIM = J
        do L=1,LLIM
          IJKL_pack = IJKL_pack+1
          if (K == L) then
            PS(IJKL_pack) = 0.5d0*PT(L,K,J,I)
          else
            PS(IJKL_pack) = 0.5d0*(PT(L,K,J,I)+PT(K,L,J,I))
            PA(IJKL_pack) = 0.5d0*(PT(L,K,J,I)-PT(K,L,J,I))
          end if
        end do
      end do
    end do
  end do
else
  ! special case for NACTEL = 1
  write(6,*) 'CheMPS2 does not allow 1 electron.'
end if

return

end subroutine CHEMPS2_DENSI_RASSCF
