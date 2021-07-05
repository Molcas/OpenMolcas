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

subroutine PCM_Cav_Hss(ToAng,nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,EF_n,EF_e,Sphere,ISphe,nOrd,Tessera,Q,Qtot,DM,HssPCM,DerMat, &
                       DerTes,DerPunt,DerRad,DerCentr,QDer,Der1,Der2,VDer)

implicit real*8(a-h,o-z)
integer IAtm(nAt)
dimension EF_n(3,*), EF_e(3,*)
dimension AtmC(3,nAt), AtmChg(nAt)
dimension Sphere(4,*), ISphe(*), nOrd(*)
dimension Tessera(4,*), Q(2,*), DM(nTs,*), Qtot(*), QDer(3,nAt,*)
dimension HssPCM(nAt3,*), DerMat(nTs,*), Der1(*), Der2(*), VDer(nTs,*)
dimension DerTes(nTs,NAt,3), DerPunt(nTs,NAt,3,3)
dimension DerRad(nS,NAt,3), DerCentr(nS,NAt,3,3)
logical DoPot, DoFld
#include "real.fh"

!***********************************************************
!
! Compute the PCM geometric contribution to gradients
!
!***********************************************************

call FZero(HssPCM,nAt3*nAt3)

! Compute the electric field on tesserae
DoPot = .false.
DoFld = .true.
call V_EF_PCM(nAt,nTs,DoPot,DoFld,AtmC,Tessera,Der2,EF_n,EF_e)

! Compute the derivative of the potential integrals
! (they should be passed to VDer_PCM below: presently
! not working, this quantity is temporarily read
! from a file in VDer_PCM).

call PotGrd(Der1,nAt3)

! Derivatives of the total electrostatic potential on tiles

call VDer_PCM(nAt,nTs,nS,AtmC,AtmChg,EF_n,EF_e,Tessera,iSphe,DerTes,DerPunt,DerRad,DerCentr,VDer)

! Derivatives of solvation charges

call DerChg(nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,DM,DerMat,Tessera,Q,Qtot,QDer,DerTes,DerPunt,DerCentr,DerRad,Der1,Der2,VDer, &
            Sphere,ISphe)
! pcm_solvent
!write(6,'(a)') 'Charge derivative wrt. atom and coord.'
!do iat=1,nat
!  do ic=1,3
!    do its=1,nts
!      write(6,'(3i5,f20.12)') iat,ic,its,QDer(ic,iat,its)
!    end do
!  end do
!end do
! pcm_solvent end

! Derivative of the cavity term

call DerCav(ToAng,nTs,nAt,nS,nAt3,Eps,Tessera,Q,Qtot,Der1,Dertes,DerPunt,DerCentr,DerRad,QDer,Sphere,iSphe,nOrd)

return

end subroutine PCM_Cav_Hss
