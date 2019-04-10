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
        use AI, only: npxAI, anAI, pAI, lb, make_parameters, lm_save
        real*8, allocatable :: x(:,:), y(:), dy(:), rl(:,:), &
                dl(:,:), mat(:,:), Iden(:,:),full_R(:,:), &
                nx(:,:),cv(:,:,:),Kv(:),pred(:),Ys(:),var(:),Rones(:), &
                sigma(:),l(:),gpred(:),hpred(:),ll(:,:)
        real*8  sb,variance,detR,lh !p
        real*8, parameter :: PI = 4.0 * atan (1.0_8), h=1e-5, eps=1e-6 ! eps avoid to become singular
        integer prev_ns,m_t,npx,counttimes
      end module globvar
