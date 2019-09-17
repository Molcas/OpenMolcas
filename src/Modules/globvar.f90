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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
      module globvar
        use AI, only: npxAI, anMd, pAI, lb, blAI, blvAI, mblAI, blaAI, &
                        blavAI
        use kriging
        real*8, allocatable ::  &
                rl(:,:,:), dl(:,:), full_Rinv(:,:), &
                full_R(:,:), nx(:,:), Kv(:), & !Iden(:,:),
                cv(:,:,:,:), cvg(:,:,:),cvh(:,:,:,:), &
                var(:), Rones(:), sigma(:), l(:), &
                pred(:), gpred(:,:), hpred(:,:,:), ll(:), &
                cvMatFder(:,:), cvMatSder(:,:), cvMatTder(:,:)
        real*8  sb,variance,detR,lh,sbO,sbmev !p
        real*8, parameter :: PI = 4.0 * atan (1.0_8), h=1e-5, &
                  eps = 1e-14, eps2 = 1e-14
! eps avoid to become singular in 1st der & eps2 in 2nd der
        integer prev_ns, m_t, npx, counttimes
      end module globvar
