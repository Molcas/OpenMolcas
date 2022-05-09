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
! Copyright (C) 1991,2021,2022, Roland Lindh                           *
!***********************************************************************

subroutine Do_NInt_d()

use nq_Grid, only: GradRho, Grid_AO, iBfn_Index, TabAO, vLapl, vRho, vSigma, vTau, Weights
use nq_Info, only: Functional_type, GGA_type, LDA_type, meta_GGA_type1, meta_GGA_type2
use Constants, only: Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iCB, iGrid, mGrid, nBfn, nD
real(kind=wp) :: gx, gxa, gxb, gy, gya, gyb, gz, gza, gzb, Temp0, Temp0a, Temp0b, Temp1, Temp1a, Temp1b, Temp2, Temp2a, Temp2b, &
                 Temp3, Temp3a, Temp3b, Temp4, Temp45, Temp45a, Temp45b, Temp4a, Temp4b, Temp5, Temp5a, Temp5b, Tmp, Tmp1, Tmp2

!                                                                      *
!***********************************************************************
!                                                                      *
nD = size(Grid_AO,4)
mGrid = size(TabAO,2)
nBfn = size(iBfn_Index,2)
!                                                                      *
!***********************************************************************
!                                                                      *
select case (Functional_type)

  case (LDA_type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    ! F(Rho)
    !
    ! for the integrals we need:
    !
    ! phi_i dF/dRho phi_j
    !
    ! Grid_AO contains
    ! 1: phi_i dF/dRho
    !
    ! Final integral assembled as, done in do_nIntx.
    !
    ! Grid_AO(1)_i phi_j
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    select case (nD)
      !                                                                *
      !*****************************************************************
      !                                                                *
      case (1)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid) < Thr) cycle
          Tmp = vRho(1,iGrid)*Weights(iGrid)

          do iCB=1,nBfn
            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Tmp
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case (2)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid)+Rho(2,iGrid) < Thr) cycle
          Tmp1 = vRho(1,iGrid)*Weights(iGrid)
          Tmp2 = vRho(2,iGrid)*Weights(iGrid)

          do iCB=1,nBfn
            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Tmp1
            Grid_AO(1,iGrid,iCB,2) = TabAO(1,iGrid,iCB)*Tmp2
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case default
        write(u6,*) 'Invalid nD value:',nD
        call Abend()
    end select  ! nD
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (GGA_type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    ! F(Rho,Sigma) : Sigma=GradRho*GradRho
    !
    ! dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
    !
    ! for the integrals we need:
    !
    !    phi_i dF/dRho phi_j
    ! +  phi_i 2 (dF/dSigma) {GradRho Grad(phi_j)}
    ! +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
    !
    ! Grid_AO contains
    !    1:  0.5 * phi_i dF/dRho
    !      + {Grad(phi_i GradRho} 2 (dF/dSigma)
    !
    ! Final integral assembled as, done in do_nIntx.
    !
    !    Grid_AO(1)_i phi_j
    !  + phi_i Grid_AO(1)_j
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    select case (nD)
      !                                                                *
      !*****************************************************************
      !                                                                *
      case (1)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid) < Thr) cycle
          gx = GradRho(1,iGrid)*Weights(iGrid)
          gy = GradRho(2,iGrid)*Weights(iGrid)
          gz = GradRho(3,iGrid)*Weights(iGrid)

          Temp0 = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp1 = gx*Two*vSigma(1,iGrid)
          Temp2 = gy*Two*vSigma(1,iGrid)
          Temp3 = gz*Two*vSigma(1,iGrid)

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0+TabAO(2,iGrid,iCB)*Temp1+TabAO(3,iGrid,iCB)*Temp2+ &
                                     TabAO(4,iGrid,iCB)*Temp3
          end do

        end do
        !                                                              *
        !**************************************************************
        !                                                              *
      case (2)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid)+Rho(2,iGrid) < Thr) cycle
          gxa = Gradrho(1,iGrid)*Weights(iGrid)
          gya = Gradrho(2,iGrid)*Weights(iGrid)
          gza = Gradrho(3,iGrid)*Weights(iGrid)
          gxb = Gradrho(4,iGrid)*Weights(iGrid)
          gyb = Gradrho(5,iGrid)*Weights(iGrid)
          gzb = Gradrho(6,iGrid)*Weights(iGrid)

          Temp0a = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp0b = Half*vRho(2,iGrid)*Weights(iGrid)
          Temp1a = Two*vSigma(1,iGrid)*gxa+vSigma(2,iGrid)*gxb
          Temp1b = Two*vSigma(3,iGrid)*gxb+vSigma(2,iGrid)*gxa
          Temp2a = Two*vSigma(1,iGrid)*gya+vSigma(2,iGrid)*gyb
          Temp2b = Two*vSigma(3,iGrid)*gyb+vSigma(2,iGrid)*gya
          Temp3a = Two*vSigma(1,iGrid)*gza+vSigma(2,iGrid)*gzb
          Temp3b = Two*vSigma(3,iGrid)*gzb+vSigma(2,iGrid)*gza

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0a+TabAO(2,iGrid,iCB)*Temp1a+TabAO(3,iGrid,iCB)*Temp2a+ &
                                     TabAO(4,iGrid,iCB)*Temp3a
            Grid_AO(1,iGrid,iCB,2) = TabAO(1,iGrid,iCB)*Temp0b+TabAO(2,iGrid,iCB)*Temp1b+TabAO(3,iGrid,iCB)*Temp2b+ &
                                     TabAO(4,iGrid,iCB)*Temp3b
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case default
        write(u6,*) 'Invalid nD value:',nD
        call Abend()
    end select  ! nD
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type1)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    ! F(Rho,Sigma,Tau) : Sigma=GradRho*GradRho
    !
    ! dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
    !
    ! for the integrals we need:
    !
    !    phi_i dF/dRho phi_j
    ! +  phi_i 2 (dF/dSigma) {GradRho Grad(phi_j)}
    ! +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
    ! +  dF/dTau {Grad(phi_i) Grad(phi_j)}
    !
    ! Grid_AO contains
    !    1: 0.5 * phi_i dF/dRho + {Grad(phi_i) GradRho} 2 (dF/dSigma)
    !    2: Grad(phi_i)_x dF/dTau
    !    3: Grad(phi_i)_y dF/dTau
    !    4: Grad(phi_i)_z dF/dTau
    !
    ! Final integral assembled as, done in do_nIntx.
    !
    !    Grid_AO(1)_i phi_j
    !  + Phi_i * Grid_AO(1)_j
    !  + Grid_AO(2)_i Grad(phi_j)_x
    !  + Grid_AO(3)_i Grad(phi_j)_y
    !  + Grid_AO(4)_i Grad(phi_j)_z
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    select case (nD)
      !                                                                *
      !*****************************************************************
      !                                                                *
      case (1)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid) < Thr) cycle
          gx = GradRho(1,iGrid)*Weights(iGrid)
          gy = GradRho(2,iGrid)*Weights(iGrid)
          gz = GradRho(3,iGrid)*Weights(iGrid)

          Temp0 = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp1 = gx*Two*vSigma(1,iGrid)
          Temp2 = gy*Two*vSigma(1,iGrid)
          Temp3 = gz*Two*vSigma(1,iGrid)

          Temp4 = Half*vTau(1,iGrid)*Weights(iGrid)

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0+TabAO(2,iGrid,iCB)*Temp1+TabAO(3,iGrid,iCB)*Temp2+ &
                                     TabAO(4,iGrid,iCB)*Temp3
            Grid_AO(2,iGrid,iCB,1) = TabAO(2,iGrid,iCB)*Temp4
            Grid_AO(3,iGrid,iCB,1) = TabAO(3,iGrid,iCB)*Temp4
            Grid_AO(4,iGrid,iCB,1) = TabAO(4,iGrid,iCB)*Temp4
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case (2)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid)+Rho(2,iGrid) < Thr) cycle
          gxa = Gradrho(1,iGrid)*Weights(iGrid)
          gya = Gradrho(2,iGrid)*Weights(iGrid)
          gza = Gradrho(3,iGrid)*Weights(iGrid)
          gxb = Gradrho(4,iGrid)*Weights(iGrid)
          gyb = Gradrho(5,iGrid)*Weights(iGrid)
          gzb = Gradrho(6,iGrid)*Weights(iGrid)

          Temp0a = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp0b = Half*vRho(2,iGrid)*Weights(iGrid)
          Temp1a = Two*vSigma(1,iGrid)*gxa+vSigma(2,iGrid)*gxb
          Temp1b = Two*vSigma(3,iGrid)*gxb+vSigma(2,iGrid)*gxa
          Temp2a = Two*vSigma(1,iGrid)*gya+vSigma(2,iGrid)*gyb
          Temp2b = Two*vSigma(3,iGrid)*gyb+vSigma(2,iGrid)*gya
          Temp3a = Two*vSigma(1,iGrid)*gza+vSigma(2,iGrid)*gzb
          Temp3b = Two*vSigma(3,iGrid)*gzb+vSigma(2,iGrid)*gza
          Temp4a = Half*vTau(1,iGrid)*Weights(iGrid)
          Temp4b = Half*vTau(2,iGrid)*Weights(iGrid)

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0a+TabAO(2,iGrid,iCB)*Temp1a+TabAO(3,iGrid,iCB)*Temp2a+ &
                                     TabAO(4,iGrid,iCB)*Temp3a
            Grid_AO(2,iGrid,iCB,1) = TabAO(2,iGrid,iCB)*Temp4a
            Grid_AO(3,iGrid,iCB,1) = TabAO(3,iGrid,iCB)*Temp4a
            Grid_AO(4,iGrid,iCB,1) = TabAO(4,iGrid,iCB)*Temp4a

            Grid_AO(1,iGrid,iCB,2) = TabAO(1,iGrid,iCB)*Temp0b+TabAO(2,iGrid,iCB)*Temp1b+TabAO(3,iGrid,iCB)*Temp2b+ &
                                     TabAO(4,iGrid,iCB)*Temp3b
            Grid_AO(2,iGrid,iCB,2) = TabAO(2,iGrid,iCB)*Temp4b
            Grid_AO(3,iGrid,iCB,2) = TabAO(3,iGrid,iCB)*Temp4b
            Grid_AO(4,iGrid,iCB,2) = TabAO(4,iGrid,iCB)*Temp4b

          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case default
        write(u6,*) 'Invalid nD value:',nD
        call Abend()
    end select  ! nD
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type2)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    ! F(Rho,Sigma,Tau,Lapl) : Sigma=GradRho*GradRho
    !
    ! dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
    !
    ! for the integrals we need:
    !
    !    phi_i dF/dRho phi_j
    ! +  phi_i 2 dF/dSigma {GradRho Grad(phi_j)}
    ! +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
    ! +  dF/dTau {Grad(phi_i) Grad(phi_j)}
    ! +  dF/dLapl Lapl(phi_i) phi_j
    ! +  dF/dLapl 2 {Grad(phi_i) Grad(phi_j)}
    ! +  dF/dLapl phi_i Lapl(phi_j)
    !
    ! Grid_AO contains
    !    1: 0.5 phi_i dF/dRho + {Grad(phi_i) GradRho} 2 (dF/dSigma)
    !      +Lapl(phi_i) dF/dLapl
    !    2: Grad(phi_i)_x (dF/dTau + 2 dF/dLapl)
    !    3: Grad(phi_i)_y (dF/dTau + 2 dF/dLapl)
    !    4: Grad(phi_i)_z (dF/dTau + 2 dF/dLapl)
    !
    ! Final integral assembled as, done in do_nIntx.
    !
    !    Grid_AO(1)_i phi_j
    !  + phi_i Grid_AO(1)_j
    !  + Grid_AO(2)_i Grad(phi_j)_x
    !  + Grid_AO(3)_i Grad(phi_j)_y
    !  + Grid_AO(4)_i Grad(phi_j)_z
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    select case (nD)
      !                                                                *
      !*****************************************************************
      !                                                                *
      case (1)
        !                                                              *
        !***************************************************************
        !                                                              *

        do iGrid=1,mGrid

          !if (Rho(1,iGrid) < Thr) cycle
          gx = GradRho(1,iGrid)*Weights(iGrid)
          gy = GradRho(2,iGrid)*Weights(iGrid)
          gz = GradRho(3,iGrid)*Weights(iGrid)

          Temp0 = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp1 = gx*Two*vSigma(1,iGrid)
          Temp2 = gy*Two*vSigma(1,iGrid)
          Temp3 = gz*Two*vSigma(1,iGrid)

          Temp4 = Half*vTau(1,iGrid)*Weights(iGrid)
          Temp5 = vLapl(1,iGrid)*Weights(iGrid)
          Temp45 = Temp4+Two*Temp5

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0+TabAO(2,iGrid,iCB)*Temp1+TabAO(3,iGrid,iCB)*Temp2+ &
                                     TabAO(4,iGrid,iCB)*Temp3+(TabAO(5,iGrid,iCB)+TabAO(8,iGrid,iCB)+TabAO(10,iGrid,iCB))*Temp5
            Grid_AO(2,iGrid,iCB,1) = TabAO(2,iGrid,iCB)*Temp45
            Grid_AO(3,iGrid,iCB,1) = TabAO(3,iGrid,iCB)*Temp45
            Grid_AO(4,iGrid,iCB,1) = TabAO(4,iGrid,iCB)*Temp45
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case (2)
        !                                                              *
        !***************************************************************
        !                                                              *
        do iGrid=1,mGrid

          !if (Rho(1,iGrid)+Rho(2,iGrid) < Thr) cycle
          gxa = Gradrho(1,iGrid)*Weights(iGrid)
          gya = Gradrho(2,iGrid)*Weights(iGrid)
          gza = Gradrho(3,iGrid)*Weights(iGrid)
          gxb = Gradrho(4,iGrid)*Weights(iGrid)
          gyb = Gradrho(5,iGrid)*Weights(iGrid)
          gzb = Gradrho(6,iGrid)*Weights(iGrid)

          Temp0a = Half*vRho(1,iGrid)*Weights(iGrid)
          Temp0b = Half*vRho(2,iGrid)*Weights(iGrid)
          Temp1a = Two*vSigma(1,iGrid)*gxa+vSigma(2,iGrid)*gxb
          Temp1b = Two*vSigma(3,iGrid)*gxb+vSigma(2,iGrid)*gxa
          Temp2a = Two*vSigma(1,iGrid)*gya+vSigma(2,iGrid)*gyb
          Temp2b = Two*vSigma(3,iGrid)*gyb+vSigma(2,iGrid)*gya
          Temp3a = Two*vSigma(1,iGrid)*gza+vSigma(2,iGrid)*gzb
          Temp3b = Two*vSigma(3,iGrid)*gzb+vSigma(2,iGrid)*gza
          Temp4a = Half*vTau(1,iGrid)*Weights(iGrid)
          Temp4b = Half*vTau(2,iGrid)*Weights(iGrid)
          Temp5a = vLapl(1,iGrid)*Weights(iGrid)
          Temp5b = vLapl(2,iGrid)*Weights(iGrid)
          Temp45a = Temp4a+Two*Temp5a
          Temp45b = Temp4b+Two*Temp5b

          do iCB=1,nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO(1,iGrid,iCB)*Temp0a+TabAO(2,iGrid,iCB)*Temp1a+TabAO(3,iGrid,iCB)*Temp2a+ &
                                     TabAO(4,iGrid,iCB)*Temp3a+(TabAO(5,iGrid,iCB)+TabAO(8,iGrid,iCB)+TabAO(10,iGrid,iCB))*Temp5a
            Grid_AO(2,iGrid,iCB,1) = TabAO(2,iGrid,iCB)*Temp45a
            Grid_AO(3,iGrid,iCB,1) = TabAO(3,iGrid,iCB)*Temp45a
            Grid_AO(4,iGrid,iCB,1) = TabAO(4,iGrid,iCB)*Temp45a

            Grid_AO(1,iGrid,iCB,2) = TabAO(1,iGrid,iCB)*Temp0b+TabAO(2,iGrid,iCB)*Temp1b+TabAO(3,iGrid,iCB)*Temp2b+ &
                                     TabAO(4,iGrid,iCB)*Temp3b+(TabAO(5,iGrid,iCB)+TabAO(8,iGrid,iCB)+TabAO(10,iGrid,iCB))*Temp5b
            Grid_AO(2,iGrid,iCB,2) = TabAO(2,iGrid,iCB)*Temp45b
            Grid_AO(3,iGrid,iCB,2) = TabAO(3,iGrid,iCB)*Temp45b
            Grid_AO(4,iGrid,iCB,2) = TabAO(4,iGrid,iCB)*Temp45b
          end do

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      case default
        write(u6,*) 'Invalid nD value:',nD
        call Abend()
    end select  ! nD
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case default
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    write(u6,*) 'DFT_Int: Illegal functional type!'
    call Abend()
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Do_NInt_d
