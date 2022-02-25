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

use OFembed, only: dFMD
use OFembed, only: Rep_EN, Func_AB, Func_A, Func_B, Energy_NAD, V_Nuc_AB, V_Nuc_BA, V_emb

implicit real*8(a-h,o-z)
#include "Molcas.fh"
real*8 ReCharge(MxAtom)
character*16 NamRfil
character*10 Fmt
integer Cho_X_GetTol
external Cho_X_GetTol

call Get_iScalar('nSym',nSym)
call Get_iScalar('Unique atoms',nAtoms)
call Get_dArray('Effective nuclear Charge',ReCharge,nAtoms)

call Get_NameRun(NamRfil)
call NameRun('AUXRFIL')
call PotNuc_nad(nSym,nAtoms,ReCharge,ZRE_nad)

call Get_dEnergy(Energy_B)
if (dFMD > 0.0) call Get_dScalar('KSDFT energy',Ec_A)

call NameRun(NamRfil)

iTol = Cho_X_GetTol(8)
call Add_Info('V_OFE',[V_emb],1,iTol)
call Add_Info('V_NUC',[V_Nuc_AB],1,iTol)
call Add_Info('E_NAD',[Energy_NAD],1,iTol)
call Add_Info('RP_EN',[Rep_EN],1,iTol)

Fmt = '(A,F19.10)'
write(6,*)
write(6,*) '     -----------------------------------------------'
write(6,*) '      Orbital-Free Embedding Calculation : Results  '
write(6,*) '     -----------------------------------------------'
write(6,Fmt) '        DFT energy  (A)    : ',Func_A
write(6,Fmt) '        DFT energy  (B)    : ',Func_B
write(6,Fmt) '        DFT energy (A+B)   : ',Func_AB
write(6,*)
write(6,Fmt) '        Nonelectr. Vemb    : ',V_emb ! for <A|Vq|A>
write(6,*)
write(6,Fmt) '        Energy (A)         : ',Energy_A
write(6,Fmt) '        Energy (B)         : ',Energy_B
write(6,Fmt) '        DFT energy (NAD)   : ',Energy_NAD
write(6,Fmt) '        Vnuc(B)*rhoA       : ',V_Nuc_AB
write(6,Fmt) '        Vnuc(A)*rhoB       : ',V_Nuc_BA
write(6,Fmt) '        Electr. repulsion  : ',Rep_EN
write(6,*) '     -----------------------------------------------'
write(6,Fmt) '       Nuclear rep. (A--B) : ',ZRE_nad
write(6,Fmt) '       Energy (A+B)        : ',Energy_B+Energy_A+Energy_NAD+V_Nuc_AB+V_Nuc_BA+Rep_EN+ZRE_nad
if (dFMD > 0.0) write(6,Fmt) '       SCF restoring Ec(A) : ',Ec_A
write(6,*) '     -----------------------------------------------'
write(6,*)
write(6,*)

call Put_dScalar('NAD dft energy',Energy_NAD)

return

end subroutine OFE_print
