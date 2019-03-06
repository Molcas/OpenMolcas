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
        real*8, allocatable :: x(:,:), y(:), dy(:), rl(:,:), dl(:,:), mat(:,:), &
                Iden(:,:),full_R(:,:), &
                nx(:),cv(:,:),lb(:),Kv(:),pred(:),Ys(:),var(:),Rones(:), &
                sigma(:),gpred(:),hpred(:)
        real*8 p,l,sb,variance,ll,detR
        real*8, parameter :: PI = 4.0 * atan (1.0_8), h=1e-5, eps=1e-6 ! eps avoid to become singular
        integer NS,dims,prev_ns,m_t,npx,counttimes
        logical :: anamat=.true. !Analitical Matern
      end module globvar