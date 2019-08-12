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
        SUBROUTINE matern(dh, m, d1, d2)
            use globvar
#include "stdalloc.fh"
            integer d1,d2,i
            REAL*8 a,d,dh(d1,d2),m(d1,d2)
            REAL*8, Allocatable :: d0(:,:)
            INTEGER*8 c
!
            Call mma_Allocate(d0,d1,d2,label="d0")
!
! For this expresion you can check https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
! and equations (11) and (12) on ref.
            d0 = sqrt(dh)
            c = idint(pAI)
            a = Gamma(pAI+1)/Gamma(2*pAI+1)
            m = 0.0D0
            do i=0, c
                d=DBLE(i)
                m = m + (Gamma(pAI+1.0D0+d)/(Gamma(d+1.D0)*Gamma(pAI+1.0D0-d)))*(2.0D0*Sqrt(2.0D0*pAI+1.0D0)*d0)**(pAI-i)
            enddo
            m = a*m*exp(-sqrt(2.0D0*pAI+1)*d0)
!
            Call mma_deallocate(d0)
!
        END

        SUBROUTINE matderiv(nd, d, m, d1, d2)
            use globvar
#include "stdalloc.fh"
            integer nd,d1,d2,p0,k
            real*8 nr,kr,a,d(d1,d2),m(d1,d2),t
            real*8, Allocatable :: b(:,:), dh(:,:), c(:,:)
!
            Call mma_Allocate(b,d1,d2,label="b")
            Call mma_Allocate(dh,d1,d2,label="dh")
            Call mma_Allocate(c,d1,d2,label="c")
!
            m = 0
            if (anMd) then
                p0=int(pAI)
                t=sqrt(2.0D0*pAI+1)
                dh = sqrt(d)
                c=(2.0D0*pAI+1)/(2.0D0*pAI-1)*exp(-t*dh)
                if (pAI.gt.3.or.pAI.lt.1) then
                    Write(6,*) 'Analytical Matern derivatives (anamat=.True.)'
                    Write(6,*) 'is only valid for pAI = 1, 2 and 3(v = 3/2, 5/2 and 7/2)'
                else
                    select case (p0)
                        case (1)
                            select case (nd)
                                case (1)
                                    m = -c/2.0D0
                                case (2)
                                    m = c*merge(0.75D0*t/dh,dh,dh.ne.0)/3.0d0
                                case (3)
                                    m = -(2.0D0*t-3.0D0*dh)*c
                            end select
                        case (2)
                            select case (nd)
                                case (1)
                                    ! m = -c*dh*(1.0+t*dh) !New way
                                    m =-c*(1.0D0+t*dh)/2.0d0
                                    ! write (6,*) '1d m',m
                                case (2)
                                    ! m = c*(dh*(5*dh-t)+1) !New way
                                    m = c*5.0D0/4.0D0
                                    ! write (6,*) '2d m',m
                                case (3)
                                    ! m = -c*dh*(5*t*dh-15) !New way
                                    m = merge(-5.0D0/8.0D0*t/dh,dh,dh.ne.0)*c
                                    ! write (6,*) '3d m',m
                            end select
                        case (3)
                            ! write (6,*) 'Analitical Matern derivatives num',nd
                            select case (nd)
                                case (1)
                                    m =-c*(1.0D0+t*dh+dh**2)/2.0D0
                                case (2)
                                    m = c*7.0D0*(1.0D0+t*dh)/12.0D0
                                case (3)
                                    m = -c*49.0D0/24.0D0
                                    ! write (6,*) '3th der dh',dh
                            end select
                    end select
                endif
            else
                ! write (6,*) 'Numerical Matern derivatives num',nd
                nr = dble(nd)
                a = Gamma(nr+1.0D0)/h**nd
                b = 0.0D0
                do k = 0, nd
                    kr = dble(k)
                    dh = d + kr*h
                    call matern(dh, m, size(dh,1), size(dh,2))
                    b = b + DBLE((-1)**(k+1))/(Gamma(nr-kr+1.0D0)*Gamma(kr+1.0D0))*m
                enddo
                m = a*b*DBLE((-1)**(nr+1))
            endif
!
            Call mma_deAllocate(b)
            Call mma_deAllocate(dh)
            Call mma_deAllocate(c)
!
        END
