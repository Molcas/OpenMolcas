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

subroutine DrvEMBg(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradient with respect to the OFE   *
!         energy (DFT-type contribution only)                          *
!                                                                      *
!                                                                      *
!         int{ [Vxc^nad(r) + Ts^nad(r) + Vnuc^B(r)] (drhoA/dRk) dr }   *
!                                                                      *
!                                                                      *
! Note:  rho=rho_A+rho_B  where the dmat elements of rho_B on the bsfs *
!                         of A have been annihilated. This also gives  *
!                         drhoA/dRk=drho/dRk (for Rk in A)             *
!                                                                      *
! Called from: Alaska                                                  *
!                                                                      *
!                                                                      *
!     Author: F. Aquilante, Dept. of Phys. Chem.                       *
!             University of Geneva, Switzerland                        *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Para_Info, only: King
use OFembed, only: OFE_KSDFT

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
character Label*80
real*8 Grad(nGrad), Temp(nGrad)
logical Do_Grad

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!...  Prologue
iRout = 131
iPrint = nPrint(iRout)
LuWr = 6

call StatusLine(' Alaska:',' Computing OFembedding gradients')

call Set_Basis_Mode('Valence')
call Setup_iSD()
nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

!     DFT-OFE    gradient                                              *
!***********************************************************************
!                                                                      *
Do_Grad = .true.
call DrvEMB_(nDens,OFE_KSDFT,Do_Grad,Temp,nGrad,'SCF ')

iEnd = 1
99 continue
if (OFE_KSDFT(iEnd:iEnd) == ' ') then
  iEnd = iEnd-1
else
  iEnd = iEnd+1
  Go To 99
end if
Label = 'DFT-OFE('//OFE_KSDFT(1:iEnd)//') contribution'
jPrint = nPrint(112)
if (jPrint >= 15) call PrGrad(Label,Temp,nGrad,ChDisp,5)
if (king()) call DaXpY_(nGrad,One,Temp,1,Grad,1)
if (iPrint < 6) Go To 777
write(LuWr,*)
if (Grid_Type == Moving_Grid) then
  write(LuWr,*) 'DFT-OFE contribution computed for a moving grid.'
else
  write(LuWr,*) 'DFT-OFE contribution computed for a fixed grid.'
end if
write(LuWr,*)
777 continue
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
call CWTime(TCpu2,TWall2)
call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)

return

end subroutine DrvEMBg
