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
SUBROUTINE kernels()
    use globvar
!    use ogpf
!    type(gpf):: gp
    integer i,z
    real*8 tpred(npx)
    character(len=100) foo
    character(len=10) foonum
    character(len=10) foovar
    character(len=10) fooana
    character(len=10) foop
!    character(len=50) footitle
    allocate (kv(m_t),pred(npx),var(npx),sigma(npx),cv(m_t,npx))
!    call gp%xlabel ('x')
!    call gp%ylabel ('y')
!    call gp%options('set label 1 at graph 0.02, 0.02 tc lt 3')
    z=int(lb(3))
    do i = 1,z
        l=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
!        print *,"New value l: ",l
        call covarmatrix()
        call k()
        call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict()
        write(foonum,'(f0.3)') l
        write(foovar,'(f0.3)') variance
        write(fooana,'(L1)') anamat
        write(foop,'(f0.1)') p
!        write(footitle,'(f0.3)') variance
!        call gp%title(footitle)
       foo='set label 1 "l: '//foonum//' Var: '//foovar//' Ana: '//fooana//' p: '//foop//'"'
!       call gp%options(foo)
!       call gp%plot(nx, pred, 't "Pred" w lines lc "blue" lt 1 lw 2','', &
!           x, y, 't "Data Points" w points pt 8','', &
!           nx, pred+sigma, 't "Pred+sigma" w lines lc "black" lt 0 lw 2','', &
!           nx, ny, 't "Orig" w lines lc "red" lt 0 lw 2')
    enddo
    tpred=pred
    call covarvector(1)
    call predict()
    gpred=pred
    call covarvector(2)
    call predict()
    hpred=pred
    pred=tpred
!    call gp%plot(nx, pred, 't "Pred" w lines lc "blue" lt 1 lw 2','', &
!        x, y, 't "Data Points" w points pt 8','', &
!        nx, gpred, 't "GEK Gradient" w lines lc "black" lt 0 lw 2')!,'', &
!        nx, hpred, 't "GEK Hessian" w lines lc "red" lt 0 lw 2')
!    call gp%plot(nx, pred, 't "Pred" w lines lc "blue" lt 1 lw 2','', &
!        x, y, 't "Data Points" w points pt 8','', &
!        nx, hpred, 't "GEK Hessian" w lines lc "black" lt 0 lw 2')!,'', &
!        nx, hpred, 't "GEK Hessian" w lines lc "red" lt 0 lw 2')

END
