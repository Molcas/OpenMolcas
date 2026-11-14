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
! Copyright (C) Giovanni Li Manni                                      *
!***********************************************************************
! History:                                                             *
! 2018 S Dong - added print outs related to scaling                    *
! 2026 Bhavnesh Jangid - Added print outs for LPDFT quadrature results *
!***********************************************************************

subroutine add_pdft_quadrature_info()

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I

implicit none

call Add_Info('dens_tt',[Dens_I],1,6)
call Add_Info('dens_a1',[Dens_a1],1,6)
call Add_Info('dens_b1',[Dens_b1],1,6)
call Add_Info('dens_a2',[Dens_a2],1,6)
call Add_Info('dens_b2',[Dens_b2],1,6)
call Add_Info('exch_f',[CoefX],1,6)
call Add_Info('corr_f',[CoefR],1,6)
call Add_Info('excha_a',[Funcaa],1,6)
call Add_Info('excha_b',[Funcbb],1,6)
call Add_Info('corr_e',[Funccc],1,6)

end subroutine add_pdft_quadrature_info
