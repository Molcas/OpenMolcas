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
        SUBROUTINE kernels(iter,nInter)
            use globvar
            integer i,z,j,iter,nInter,lm
!
            call miden(iter)
            z=int(lb(3))
!
!To be change for the optmization of the l's (the right width of the Mat'ern function)
            do i = 1,z
!In this particullary case the l(j) it does not depend on the dimensionality
!It is the same.
!Could be the case that depends on the dimesionality then for every dimension the
!l changes
                do j = 1,nInter
                    l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                enddo
!
                call covarmatrix(iter,nInter)
                call k(iter)
                ll(i)=lh
            enddo
!
            lm = MaxLoc(ll,dim=nInter)
            do j = 1,nInter
                l(j)=lb(1)+(lm-1)*(lb(2)-lb(1))/(lb(3)-1)
            enddo
            Call covarmatrix(iter,nInter)
            Call k(iter)
            Write(6,*) 'Likelihood function pred index lm,l,ll(lm),variance',lm,ll(lm),l,variance
        END

        subroutine miden(iter)
            use globvar
            integer j,iter
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine
