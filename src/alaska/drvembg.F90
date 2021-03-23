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
      SubRoutine DrvEMBg(Grad,Temp,nGrad)
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
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
      Character Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
      Logical Do_Grad
!                                                                      *
!***********************************************************************
!                                                                      *
      Call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      LuWr=6
!
      Call StatusLine(' Alaska:',' Computing OFembedding gradients')
!
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
!
!     DFT-OFE    gradient                                              *
!***********************************************************************
!                                                                      *
      Do_Grad=.True.
      Call DrvEMB_(nDens,OFE_KSDFT,Do_Grad,Temp,nGrad,'SCF ')
!
      iEnd=1
 99   Continue
      If (OFE_KSDFT(iEnd:iEnd).eq.' ') Then
         iEnd = iEnd - 1
      Else
         iEnd = iEnd + 1
         Go To 99
      End If
      Label='DFT-OFE('//OFE_KSDFT(1:iEnd)//') contribution'
      jPrint=nPrint(112)
      If (jPrint.ge.15) Call PrGrad(Label,Temp,nGrad,ChDisp,5)
      If (king()) Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      If (iPrint.lt.6) Go To 777
      Write (LuWr,*)
      If (Grid_Type.eq.Moving_Grid) Then
         Write(LuWr,*)'DFT-OFE contribution computed for a moving grid.'
      Else
         Write(LuWr,*)'DFT-OFE contribution computed for a fixed grid.'
      End If
      Write (LuWr,*)
 777  Continue
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Free_iSD()
      Call CWTime(TCpu2,TWall2)
      Call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
      Return
      End
