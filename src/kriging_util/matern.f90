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
        SUBROUTINE matern(dh,d1,d2)
            use globvar
            integer d1,d2
            REAL*8 b(d1,d2),a,d,d0(d1,d2),dh(d1,d2)
            INTEGER*8 c,i
            d0=sqrt(dh)
            c=int(p)
            a=Gamma(p+1)/Gamma(2*p+1)
            b=0.0
            do i=0, c
                d=real(i)
                b = b + (Gamma(p+1.0+d)/(Gamma(d+1)*Gamma(p+1.0-d)))*(2.0*Sqrt(2.0*p+1.0)*d0)**(p-i)
            enddo
                mat=a*b*exp(-sqrt(2*p+1)*d0)
        END

        SUBROUTINE matderiv(nd,d1,d2)
            use globvar
            integer nd,d1,d2,p0
            real*8 nr,kr,a,b(d1,d2),dh(d1,d2),c(d1,d2),t
            nr=real(nd)
            if (anamat) then
                p0=int(p)
                t=sqrt(2.0*p+1)
                dh=sqrt(dl)
                c=(2.0*p+1)/(2.0*p-1)*exp(-t*dh)
                if (p.gt.2.or.p.lt.1) then
                    Write(6,*) 'Analytical Matern derivatives (anamat=.True.) is only valid for p = 1 or 2 (v = 3/2 or 5/2)'
                else
                    select case (p0)
                        case (1)
                            select case (nd)
                                case (1)
                                    mat=-c/2.0
                                case (2)
                                    mat=c*merge(0.75*t/dh,dh,dh.ne.0)/3.0
                                case (3)
                                    mat=-(2.0*t-3.0*rl)*c
                            end select
                        case (2)
                            select case (nd)
                                case (1)
                                    mat=-c*(1.0+t*dh)/2.0
                                case (2)
                                    mat=c*5.0/4.0
                                case (3)
                                    mat=merge(-5.0/8.0*t/dh,dh,dh.ne.0)*c
                            end select
                    end select
                endif
            else
                a=Gamma(nr+1.0)/h**nd
                b=0.0
                do k=0,nd
                    kr=real(k)
                    dh=dl+kr*h
                    call matern(dh,size(dh,1),size(dh,2))
                    b=b+(-1)**(k+1)/(Gamma(nr-kr+1.0)*Gamma(kr+1.0))*mat
                enddo
                mat=a*b*(-1)**(nr+1)
            endif
        END
