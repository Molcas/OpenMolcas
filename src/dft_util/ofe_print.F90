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

subroutine OFE_print(Energy_A)

use OFembed, only: dFMD, Energy_NAD, Func_A, Func_AB, Func_B, Rep_EN, V_emb, V_Nuc_AB, V_Nuc_BA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Energy_A
integer(kind=iwp) :: iTol, nAtoms, nSym
real(kind=wp) :: Ec_A, Energy_B, ZRE_nad
real(kind=wp), allocatable :: ReCharge(:)
integer(kind=iwp), external :: Cho_X_GetTol

call Get_iScalar('nSym',nSym)
call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(ReCharge,nAtoms,label='ReCharge')
call Get_dArray('Effective nuclear Charge',ReCharge,nAtoms)

call NameRun('AUXRFIL')
call PotNuc_nad(nSym,nAtoms,ReCharge,ZRE_nad)
call mma_deallocate(ReCharge)

call Get_dEnergy(Energy_B)
if (dFMD > Zero) call Get_dScalar('KSDFT energy',Ec_A)

call NameRun('#Pop')

iTol = Cho_X_GetTol(8)
call Add_Info('V_OFE',[V_emb],1,iTol)
call Add_Info('V_NUC',[V_Nuc_AB],1,iTol)
call Add_Info('E_NAD',[Energy_NAD],1,iTol)
call Add_Info('RP_EN',[Rep_EN],1,iTol)

write(u6,*)
write(u6,*) '     -----------------------------------------------'
write(u6,*) '      Orbital-Free Embedding Calculation : Results  '
write(u6,*) '     -----------------------------------------------'
write(u6,100) '        DFT energy  (A)    : ',Func_A
write(u6,100) '        DFT energy  (B)    : ',Func_B
write(u6,100) '        DFT energy (A+B)   : ',Func_AB
write(u6,*)
write(u6,100) '        Nonelectr. Vemb    : ',V_emb ! for <A|Vq|A>
write(u6,*)
write(u6,100) '        Energy (A)         : ',Energy_A
write(u6,100) '        Energy (B)         : ',Energy_B
write(u6,100) '        DFT energy (NAD)   : ',Energy_NAD
write(u6,100) '        Vnuc(B)*rhoA       : ',V_Nuc_AB
write(u6,100) '        Vnuc(A)*rhoB       : ',V_Nuc_BA
write(u6,100) '        Electr. repulsion  : ',Rep_EN
write(u6,*) '     -----------------------------------------------'
write(u6,100) '       Nuclear rep. (A--B) : ',ZRE_nad
write(u6,100) '       Energy (A+B)        : ',Energy_B+Energy_A+Energy_NAD+V_Nuc_AB+V_Nuc_BA+Rep_EN+ZRE_nad
if (dFMD > Zero) write(u6,100) '       SCF restoring Ec(A) : ',Ec_A
write(u6,*) '     -----------------------------------------------'
write(u6,*)
write(u6,*)

call Put_dScalar('NAD dft energy',Energy_NAD)

return

100 format(A,F19.10)

end subroutine OFE_print
