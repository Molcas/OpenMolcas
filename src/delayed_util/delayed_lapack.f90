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

subroutine dbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, ldc, ldu, ldvt, n, ncc, ncvt, nru
  real*8 :: c( ldc, * ), d( * ), e( * ), u( ldu, * ), vt( ldvt, * ), work( * )
  call lb_dbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info )
end subroutine dbdsqr

subroutine dgebak( job, side, n, ilo, ihi, scale, m, v, ldv, info )
  use link_blas
  implicit none
  character :: job, side
  integer :: ihi, ilo, info, ldv, m, n
  real*8 :: scale( * ), v( ldv, * )
  call lb_dgebak( job, side, n, ilo, ihi, scale, m, v, ldv, info )
end subroutine dgebak

subroutine dgebal( job, n, a, lda, ilo, ihi, scale, info )
  use link_blas
  implicit none
  character :: job
  integer :: ihi, ilo, info, lda, n
  real*8 :: a( lda, * ), scale( * )
  call lb_dgebal( job, n, a, lda, ilo, ihi, scale, info )
end subroutine dgebal

subroutine dgebd2( m, n, a, lda, d, e, tauq, taup, work, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  real*8 :: a( lda, * ), d( * ), e( * ), taup( * ), tauq( * ), work( * )
  call lb_dgebd2( m, n, a, lda, d, e, tauq, taup, work, info )
end subroutine dgebd2

subroutine dgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, lda, lwork, m, n
  real*8 :: a( lda, * ), d( * ), e( * ), taup( * ), tauq( * ), work( * )
  call lb_dgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info )
end subroutine dgebrd

subroutine dgecon( norm, n, a, lda, anorm, rcond, work, iwork, info )
  use link_blas
  implicit none
  character :: norm
  integer :: info, lda, n
  real*8 :: anorm, rcond
  integer :: iwork( * )
  real*8 :: a( lda, * ), work( * )
  call lb_dgecon( norm, n, a, lda, anorm, rcond, work, iwork, info )
end subroutine dgecon

subroutine dgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info )
  use link_blas
  implicit none
  character :: jobvl, jobvr
  integer :: info, lda, ldvl, ldvr, lwork, n
  real*8 :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), wi( * ), work( * ), wr( * )
  call lb_dgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info )
end subroutine dgeev

subroutine dgehd2( n, ilo, ihi, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: ihi, ilo, info, lda, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgehd2( n, ilo, ihi, a, lda, tau, work, info )
end subroutine dgehd2

subroutine dgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: ihi, ilo, info, lda, lwork, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info )
end subroutine dgehrd

subroutine dgelq2( m, n, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgelq2( m, n, a, lda, tau, work, info )
end subroutine dgelq2

subroutine dgelqf( m, n, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgelqf( m, n, a, lda, tau, work, lwork, info )
end subroutine dgelqf

subroutine dgels( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
  use link_blas
  implicit none
  character :: trans
  integer :: info, lda, ldb, lwork, m, n, nrhs
  real*8 :: a( lda, * ), b( ldb, * ), work( * )
  call lb_dgels( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
end subroutine dgels

subroutine dgeqr2( m, n, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgeqr2( m, n, a, lda, tau, work, info )
end subroutine dgeqr2

subroutine dgeqrf( m, n, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dgeqrf( m, n, a, lda, tau, work, lwork, info )
end subroutine dgeqrf

subroutine dgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info )
  use link_blas
  implicit none
  character :: jobu, jobvt
  integer :: info, lda, ldu, ldvt, lwork, m, n
  real*8 :: a( lda, * ), s( * ), u( ldu, * ), vt( ldvt, * ), work( * )
  call lb_dgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info )
end subroutine dgesvd

subroutine dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
  use link_blas
  implicit none
  integer :: info, lda, ldb, n, nrhs
  integer :: ipiv( * )
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
end subroutine dgesv

subroutine dgetrf( m, n, a, lda, ipiv, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  integer :: ipiv( * )
  real*8 :: a( lda, * )
  call lb_dgetrf( m, n, a, lda, ipiv, info )
end subroutine dgetrf

recursive subroutine dgetrf2( m, n, a, lda, ipiv, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  integer :: ipiv( * )
  real*8 :: a( lda, * )
  call lb_dgetrf2( m, n, a, lda, ipiv, info )
end subroutine dgetrf2

subroutine dgetri( n, a, lda, ipiv, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, lda, lwork, n
  integer :: ipiv( * )
  real*8 :: a( lda, * ), work( * )
  call lb_dgetri( n, a, lda, ipiv, work, lwork, info )
end subroutine dgetri

subroutine dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
  use link_blas
  implicit none
  character :: trans
  integer :: info, lda, ldb, n, nrhs
  integer :: ipiv( * )
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
end subroutine dgetrs

subroutine dhseqr( job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info )
  use link_blas
  implicit none
  integer :: ihi, ilo, info, ldh, ldz, lwork, n
  character :: compz, job
  real*8 :: h( ldh, * ), wi( * ), work( * ), wr( * ), z( ldz, * )
  call lb_dhseqr( job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info )
end subroutine dhseqr

function disnan( din )
  use link_blas
  implicit none
  real*8 :: din
  logical :: disnan
  disnan=lb_disnan( din )
end function disnan

subroutine dlabad( small, large )
  use link_blas
  implicit none
  real*8 :: large, small
  call lb_dlabad( small, large )
end subroutine dlabad

subroutine dlabrd( m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy )
  use link_blas
  implicit none
  integer :: lda, ldx, ldy, m, n, nb
  real*8 :: a( lda, * ), d( * ), e( * ), taup( * ), tauq( * ), x( ldx, * ), y( ldy, * )
  call lb_dlabrd( m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy )
end subroutine dlabrd

subroutine dlacn2( n, v, x, isgn, est, kase, isave )
  use link_blas
  implicit none
  integer :: kase, n
  real*8 :: est
  integer :: isgn( * ), isave( 3 )
  real*8 :: v( * ), x( * )
  call lb_dlacn2( n, v, x, isgn, est, kase, isave )
end subroutine dlacn2

subroutine dlacpy( uplo, m, n, a, lda, b, ldb )
  use link_blas
  implicit none
  character :: uplo
  integer :: lda, ldb, m, n
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dlacpy( uplo, m, n, a, lda, b, ldb )
end subroutine dlacpy

subroutine dladiv( a, b, c, d, p, q )
  use link_blas
  implicit none
  real*8 :: a, b, c, d, p, q
  call lb_dladiv( a, b, c, d, p, q )
end subroutine dladiv

function dladiv2( a, b, c, d, r, t )
  use link_blas
  implicit none
  real*8 :: a, b, c, d, r, t
  real*8 :: dladiv2
  dladiv2=lb_dladiv2( a, b, c, d, r, t )
end function dladiv2

subroutine dladiv1( a, b, c, d, p, q )
  use link_blas
  implicit none
  real*8 :: a, b, c, d, p, q
  call lb_dladiv1( a, b, c, d, p, q )
end subroutine dladiv1

subroutine dlae2( a, b, c, rt1, rt2 )
  use link_blas
  implicit none
  real*8 :: a, b, c, rt1, rt2
  call lb_dlae2( a, b, c, rt1, rt2 )
end subroutine dlae2

subroutine dlaebz( ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info )
  use link_blas
  implicit none
  integer :: ijob, info, minp, mmax, mout, n, nbmin, nitmax
  real*8 :: abstol, pivmin, reltol
  integer :: iwork( * ), nab( mmax, * ), nval( * )
  real*8 :: ab( mmax, * ), c( * ), d( * ), e( * ), e2( * ), work( * )
  call lb_dlaebz( ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info )
end subroutine dlaebz

subroutine dlaev2( a, b, c, rt1, rt2, cs1, sn1 )
  use link_blas
  implicit none
  real*8 :: a, b, c, cs1, rt1, rt2, sn1
  call lb_dlaev2( a, b, c, rt1, rt2, cs1, sn1 )
end subroutine dlaev2

subroutine dlaexc( wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info )
  use link_blas
  implicit none
  logical :: wantq
  integer :: info, j1, ldq, ldt, n, n1, n2
  real*8 :: q( ldq, * ), t( ldt, * ), work( * )
  call lb_dlaexc( wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info )
end subroutine dlaexc

subroutine dlagtf( n, a, lambda, b, c, tol, d, in, info )
  use link_blas
  implicit none
  integer :: info, n
  real*8 :: lambda, tol
  integer :: in( * )
  real*8 :: a( * ), b( * ), c( * ), d( * )
  call lb_dlagtf( n, a, lambda, b, c, tol, d, in, info )
end subroutine dlagtf

subroutine dlagts( job, n, a, b, c, d, in, y, tol, info )
  use link_blas
  implicit none
  integer :: info, job, n
  real*8 :: tol
  integer :: in( * )
  real*8 :: a( * ), b( * ), c( * ), d( * ), y( * )
  call lb_dlagts( job, n, a, b, c, d, in, y, tol, info )
end subroutine dlagts

subroutine dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info )
  use link_blas
  implicit none
  integer :: ihi, ihiz, ilo, iloz, info, ldh, ldz, n
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), wi( * ), wr( * ), z( ldz, * )
  call lb_dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info )
end subroutine dlahqr

subroutine dlahr2( n, k, nb, a, lda, tau, t, ldt, y, ldy )
  use link_blas
  implicit none
  integer :: k, lda, ldt, ldy, n, nb
  real*8 :: a( lda, * ), t( ldt, nb ), tau( nb ), y( ldy, nb )
  call lb_dlahr2( n, k, nb, a, lda, tau, t, ldt, y, ldy )
end subroutine dlahr2

function dlaisnan( din1, din2 )
  use link_blas
  implicit none
  real*8 :: din1, din2
  logical :: dlaisnan
  dlaisnan=lb_dlaisnan( din1, din2 )
end function dlaisnan

subroutine dlaln2( ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info )
  use link_blas
  implicit none
  logical :: ltrans
  integer :: info, lda, ldb, ldx, na, nw
  real*8 :: ca, d1, d2, scale, smin, wi, wr, xnorm
  real*8 :: a( lda, * ), b( ldb, * ), x( ldx, * )
  call lb_dlaln2( ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info )
end subroutine dlaln2

function dlamch( cmach )
  use link_blas
  implicit none
  character :: cmach
  real*8 :: dlamch
  dlamch=lb_dlamch( cmach )
end function dlamch

function dlaneg( n, d, lld, sigma, pivmin, r )
  use link_blas
  implicit none
  integer :: n, r
  real*8 :: pivmin, sigma
  real*8 :: d( * ), lld( * )
  integer :: dlaneg
  dlaneg=lb_dlaneg( n, d, lld, sigma, pivmin, r )
end function dlaneg

function dlange( norm, m, n, a, lda, work )
  use link_blas
  implicit none
  character :: norm
  integer :: lda, m, n
  real*8 :: a( lda, * ), work( * )
  real*8 :: dlange
  dlange=lb_dlange( norm, m, n, a, lda, work )
end function dlange

function dlansp( norm, uplo, n, ap, work )
  use link_blas
  implicit none
  character :: norm, uplo
  integer :: n
  real*8 :: ap( * ), work( * )
  real*8 :: dlansp
  dlansp=lb_dlansp( norm, uplo, n, ap, work )
end function dlansp

function dlanst( norm, n, d, e )
  use link_blas
  implicit none
  character :: norm
  integer :: n
  real*8 :: d( * ), e( * )
  real*8 :: dlanst
  dlanst=lb_dlanst( norm, n, d, e )
end function dlanst

function dlansy( norm, uplo, n, a, lda, work )
  use link_blas
  implicit none
  character :: norm, uplo
  integer :: lda, n
  real*8 :: a( lda, * ), work( * )
  real*8 :: dlansy
  dlansy=lb_dlansy( norm, uplo, n, a, lda, work )
end function dlansy

subroutine dlanv2( a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn )
  use link_blas
  implicit none
  real*8 :: a, b, c, cs, d, rt1i, rt1r, rt2i, rt2r, sn
  call lb_dlanv2( a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn )
end subroutine dlanv2

function dlapy2( x, y )
  use link_blas
  implicit none
  real*8 :: x, y
  real*8 :: dlapy2
  dlapy2=lb_dlapy2( x, y )
end function dlapy2

function dlapy3( x, y, z )
  use link_blas
  implicit none
  real*8 :: x, y, z
  real*8 :: dlapy3
  dlapy3=lb_dlapy3( x, y, z )
end function dlapy3

subroutine dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info )
  use link_blas
  implicit none
  integer :: ihi, ihiz, ilo, iloz, info, ldh, ldz, lwork, n
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), wi( * ), work( * ), wr( * ), z( ldz, * )
  call lb_dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info )
end subroutine dlaqr0

subroutine dlaqr1( n, h, ldh, sr1, si1, sr2, si2, v )
  use link_blas
  implicit none
  real*8 :: si1, si2, sr1, sr2
  integer :: ldh, n
  real*8 :: h( ldh, * ), v( * )
  call lb_dlaqr1( n, h, ldh, sr1, si1, sr2, si2, v )
end subroutine dlaqr1

subroutine dlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, &
                   work, lwork )
  use link_blas
  implicit none
  integer :: ihiz, iloz, kbot, ktop, ldh, ldt, ldv, ldwv, ldz, lwork, n, nd, nh, ns, nv, nw
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), si( * ), sr( * ), t( ldt, * ), v( ldv, * ), work( * ), wv( ldwv, * ), z( ldz, * )
  call lb_dlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, &
                  work, lwork )
end subroutine dlaqr2

subroutine dlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, &
                   work, lwork )
  use link_blas
  implicit none
  integer :: ihiz, iloz, kbot, ktop, ldh, ldt, ldv, ldwv, ldz, lwork, n, nd, nh, ns, nv, nw
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), si( * ), sr( * ), t( ldt, * ), v( ldv, * ), work( * ), wv( ldwv, * ), z( ldz, * )
  call lb_dlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, &
                  work, lwork )
end subroutine dlaqr3

subroutine dlaqr4( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info )
  use link_blas
  implicit none
  integer :: ihi, ihiz, ilo, iloz, info, ldh, ldz, lwork, n
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), wi( * ), work( * ), wr( * ), z( ldz, * )
  call lb_dlaqr4( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info )
end subroutine dlaqr4

subroutine dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, &
                   nh, wh, ldwh )
  use link_blas
  implicit none
  integer :: ihiz, iloz, kacc22, kbot, ktop, ldh, ldu, ldv, ldwh, ldwv, ldz, n, nh, nshfts, nv
  logical :: wantt, wantz
  real*8 :: h( ldh, * ), si( * ), sr( * ), u( ldu, * ), v( ldv, * ), wh( ldwh, * ), wv( ldwv, * ), z( ldz, * )
  call lb_dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, &
                  nh, wh, ldwh )
end subroutine dlaqr5

subroutine dlar1v( n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, &
                   rqcorr, work )
  use link_blas
  implicit none
  logical :: wantnc
  integer :: b1, bn, n, negcnt, r
  real*8 :: gaptol, lambda, mingma, nrminv, pivmin, resid, rqcorr, ztz
  integer :: isuppz( * )
  real*8 :: d( * ), l( * ), ld( * ), lld( * ), work( * )
  real*8 :: z( * )
  call lb_dlar1v( n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, &
                  rqcorr, work )
end subroutine dlar1v

subroutine dlarfb( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork )
  use link_blas
  implicit none
  character :: direct, side, storev, trans
  integer :: k, ldc, ldt, ldv, ldwork, m, n
  real*8 :: c( ldc, * ), t( ldt, * ), v( ldv, * ), work( ldwork, * )
  call lb_dlarfb( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork )
end subroutine dlarfb

subroutine dlarf( side, m, n, v, incv, tau, c, ldc, work )
  use link_blas
  implicit none
  character :: side
  integer :: incv, ldc, m, n
  real*8 :: tau
  real*8 :: c( ldc, * ), v( * ), work( * )
  call lb_dlarf( side, m, n, v, incv, tau, c, ldc, work )
end subroutine dlarf

subroutine dlarfg( n, alpha, x, incx, tau )
  use link_blas
  implicit none
  integer :: incx, n
  real*8 :: alpha, tau
  real*8 :: x( * )
  call lb_dlarfg( n, alpha, x, incx, tau )
end subroutine dlarfg

subroutine dlarft( direct, storev, n, k, v, ldv, tau, t, ldt )
  use link_blas
  implicit none
  character :: direct, storev
  integer :: k, ldt, ldv, n
  real*8 :: t( ldt, * ), tau( * ), v( ldv, * )
  call lb_dlarft( direct, storev, n, k, v, ldv, tau, t, ldt )
end subroutine dlarft

subroutine dlarfx( side, m, n, v, tau, c, ldc, work )
  use link_blas
  implicit none
  character :: side
  integer :: ldc, m, n
  real*8 :: tau
  real*8 :: c( ldc, * ), v( * ), work( * )
  call lb_dlarfx( side, m, n, v, tau, c, ldc, work )
end subroutine dlarfx

subroutine dlarnv( idist, iseed, n, x )
  use link_blas
  implicit none
  integer :: idist, n
  integer :: iseed( 4 )
  real*8 :: x( * )
  call lb_dlarnv( idist, iseed, n, x )
end subroutine dlarnv

subroutine dlarra( n, d, e, e2, spltol, tnrm, nsplit, isplit, info )
  use link_blas
  implicit none
  integer :: info, n, nsplit
  real*8 :: spltol, tnrm
  integer :: isplit( * )
  real*8 :: d( * ), e( * ), e2( * )
  call lb_dlarra( n, d, e, e2, spltol, tnrm, nsplit, isplit, info )
end subroutine dlarra

subroutine dlarrb( n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info )
  use link_blas
  implicit none
  integer :: ifirst, ilast, info, n, offset, twist
  real*8 :: pivmin, rtol1, rtol2, spdiam
  integer :: iwork( * )
  real*8 :: d( * ), lld( * ), w( * ), werr( * ), wgap( * ), work( * )
  call lb_dlarrb( n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info )
end subroutine dlarrb

subroutine dlarrc( jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info )
  use link_blas
  implicit none
  character :: jobt
  integer :: eigcnt, info, lcnt, n, rcnt
  real*8 :: pivmin, vl, vu
  real*8 :: d( * ), e( * )
  call lb_dlarrc( jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info )
end subroutine dlarrc

subroutine dlarrd( range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, &
                   indexw, work, iwork, info )
  use link_blas
  implicit none
  character :: order, range
  integer :: il, info, iu, m, n, nsplit
  real*8 :: pivmin, reltol, vl, vu, wl, wu
  integer :: iblock( * ), indexw( * ), isplit( * ), iwork( * )
  real*8 :: d( * ), e( * ), e2( * ), gers( * ), w( * ), werr( * ), work( * )
  call lb_dlarrd( range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, &
                  indexw, work, iwork, info )
end subroutine dlarrd

subroutine dlarre( range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, &
                   gers, pivmin, work, iwork, info )
  use link_blas
  implicit none
  character :: range
  integer :: il, info, iu, m, n, nsplit
  real*8 :: pivmin, rtol1, rtol2, spltol, vl, vu
  integer :: iblock( * ), isplit( * ), iwork( * ), indexw( * )
  real*8 :: d( * ), e( * ), e2( * ), gers( * ), w( * ),werr( * ), wgap( * ), work( * )
  call lb_dlarre( range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, &
                  gers, pivmin, work, iwork, info )
end subroutine dlarre

subroutine dlarrf( n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info )
  use link_blas
  implicit none
  integer :: clstrt, clend, info, n
  real*8 :: clgapl, clgapr, pivmin, sigma, spdiam
  real*8 :: d( * ), dplus( * ), l( * ), ld( * ), lplus( * ), w( * ), wgap( * ), werr( * ), work( * )
  call lb_dlarrf( n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info )
end subroutine dlarrf

subroutine dlarrj( n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info )
  use link_blas
  implicit none
  integer :: ifirst, ilast, info, n, offset
  real*8 :: pivmin, rtol, spdiam
  integer :: iwork( * )
  real*8 :: d( * ), e2( * ), w( * ), werr( * ), work( * )
  call lb_dlarrj( n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info )
end subroutine dlarrj

subroutine dlarrk( n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
  use link_blas
  implicit none
  integer :: info, iw, n
  real*8 :: pivmin, reltol, gl, gu, w, werr
  real*8 :: d( * ), e2( * )
  call lb_dlarrk( n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
end subroutine dlarrk

subroutine dlarrr( n, d, e, info )
  use link_blas
  implicit none
  integer :: n, info
  real*8 :: d( * ), e( * )
  call lb_dlarrr( n, d, e, info )
end subroutine dlarrr

subroutine dlarrv( n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, &
                   ldz, isuppz, work, iwork, info )
  use link_blas
  implicit none
  integer :: dol, dou, info, ldz, m, n
  real*8 :: minrgp, pivmin, rtol1, rtol2, vl, vu
  integer :: iblock( * ), indexw( * ), isplit( * ), isuppz( * ), iwork( * )
  real*8 :: d( * ), gers( * ), l( * ), w( * ), werr( * ), wgap( * ), work( * )
  real*8 :: z( ldz, * )
  call lb_dlarrv( n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, &
                  isuppz, work, iwork, info )
end subroutine dlarrv

subroutine dlartg( f, g, cs, sn, r )
  use link_blas
  implicit none
  real*8 :: cs, f, g, r, sn
  call lb_dlartg( f, g, cs, sn, r )
end subroutine dlartg

subroutine dlaruv( iseed, n, x )
  use link_blas
  implicit none
  integer :: n
  integer :: iseed( 4 )
  real*8 :: x( n )
  call lb_dlaruv( iseed, n, x )
end subroutine dlaruv

subroutine dlas2( f, g, h, ssmin, ssmax )
  use link_blas
  implicit none
  real*8 :: f, g, h, ssmax, ssmin
  call lb_dlas2( f, g, h, ssmin, ssmax )
end subroutine dlas2

subroutine dlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
  use link_blas
  implicit none
  character :: type
  integer :: info, kl, ku, lda, m, n
  real*8 :: cfrom, cto
  real*8 :: a( lda, * )
  call lb_dlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
end subroutine dlascl

subroutine dlaset( uplo, m, n, alpha, beta, a, lda )
  use link_blas
  implicit none
  character :: uplo
  integer :: lda, m, n
  real*8 :: alpha, beta
  real*8 :: a( lda, * )
  call lb_dlaset( uplo, m, n, alpha, beta, a, lda )
end subroutine dlaset

subroutine dlasq1( n, d, e, work, info )
  use link_blas
  implicit none
  integer :: info, n
  real*8 :: d( * ), e( * ), work( * )
  call lb_dlasq1( n, d, e, work, info )
end subroutine dlasq1

subroutine dlasq2( n, z, info )
  use link_blas
  implicit none
  integer :: info, n
  real*8 :: z( * )
  call lb_dlasq2( n, z, info )
end subroutine dlasq2

subroutine dlasq3( i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau )
  use link_blas
  implicit none
  logical :: ieee
  integer :: i0, iter, n0, ndiv, nfail, pp, ttype
  real*8 :: desig, dmin, dmin1, dmin2, dn, dn1, dn2, g, qmax, sigma, tau
  real*8 :: z( * )
  call lb_dlasq3( i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau )
end subroutine dlasq3

subroutine dlasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g )
  use link_blas
  implicit none
  integer :: i0, n0, n0in, pp, ttype
  real*8 :: dmin, dmin1, dmin2, dn, dn1, dn2, g, tau
  real*8 :: z( * )
  call lb_dlasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g )
end subroutine dlasq4

subroutine dlasq5( i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps )
  use link_blas
  implicit none
  logical :: ieee
  integer :: i0, n0, pp
  real*8 :: dmin, dmin1, dmin2, dn, dnm1, dnm2, tau, sigma, eps
  real*8 :: z( * )
  call lb_dlasq5( i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps )
end subroutine dlasq5

subroutine dlasq6( i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2 )
  use link_blas
  implicit none
  integer :: i0, n0, pp
  real*8 :: dmin, dmin1, dmin2, dn, dnm1, dnm2
  real*8 :: z( * )
  call lb_dlasq6( i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2 )
end subroutine dlasq6

subroutine dlasr( side, pivot, direct, m, n, c, s, a, lda )
  use link_blas
  implicit none
  character :: direct, pivot, side
  integer :: lda, m, n
  real*8 :: a( lda, * ), c( * ), s( * )
  call lb_dlasr( side, pivot, direct, m, n, c, s, a, lda )
end subroutine dlasr

subroutine dlasrt( id, n, d, info )
  use link_blas
  implicit none
  character :: id
  integer :: info, n
  real*8 :: d( * )
  call lb_dlasrt( id, n, d, info )
end subroutine dlasrt

subroutine dlassq( n, x, incx, scale, sumsq )
  use link_blas
  implicit none
  integer :: incx, n
  real*8 :: scale, sumsq
  real*8 :: x( * )
  call lb_dlassq( n, x, incx, scale, sumsq )
end subroutine dlassq

subroutine dlasv2( f, g, h, ssmin, ssmax, snr, csr, snl, csl )
  use link_blas
  implicit none
  real*8 :: csl, csr, f, g, h, snl, snr, ssmax, ssmin
  call lb_dlasv2( f, g, h, ssmin, ssmax, snr, csr, snl, csl )
end subroutine dlasv2

subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
  use link_blas
  implicit none
  integer :: incx, k1, k2, lda, n
  integer :: ipiv( * )
  real*8 :: a( lda, * )
  call lb_dlaswp( n, a, lda, k1, k2, ipiv, incx )
end subroutine dlaswp

subroutine dlasy2( ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info )
  use link_blas
  implicit none
  logical :: ltranl, ltranr
  integer :: info, isgn, ldb, ldtl, ldtr, ldx, n1, n2
  real*8 :: scale, xnorm
  real*8 :: b( ldb, * ), tl( ldtl, * ), tr( ldtr, * ), x( ldx, * )
  call lb_dlasy2( ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info )
end subroutine dlasy2

subroutine dlatrd( uplo, n, nb, a, lda, e, tau, w, ldw )
  use link_blas
  implicit none
  character :: uplo
  integer :: lda, ldw, n, nb
  real*8 :: a( lda, * ), e( * ), tau( * ), w( ldw, * )
  call lb_dlatrd( uplo, n, nb, a, lda, e, tau, w, ldw )
end subroutine dlatrd

subroutine dlatrs( uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info )
  use link_blas
  implicit none
  character :: diag, normin, trans, uplo
  integer :: info, lda, n
  real*8 :: scale
  real*8 :: a( lda, * ), cnorm( * ), x( * )
  call lb_dlatrs( uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info )
end subroutine dlatrs

subroutine dopgtr( uplo, n, ap, tau, q, ldq, work, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, ldq, n
  real*8 :: ap( * ), q( ldq, * ), tau( * ), work( * )
  call lb_dopgtr( uplo, n, ap, tau, q, ldq, work, info )
end subroutine dopgtr

subroutine dopmtr( side, uplo, trans, m, n, ap, tau, c, ldc, work, info )
  use link_blas
  implicit none
  character :: side, trans, uplo
  integer :: info, ldc, m, n
  real*8 :: ap( * ), c( ldc, * ), tau( * ), work( * )
  call lb_dopmtr( side, uplo, trans, m, n, ap, tau, c, ldc, work, info )
end subroutine dopmtr

subroutine dorg2l( m, n, k, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, k, lda, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorg2l( m, n, k, a, lda, tau, work, info )
end subroutine dorg2l

subroutine dorg2r( m, n, k, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, k, lda, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorg2r( m, n, k, a, lda, tau, work, info )
end subroutine dorg2r

subroutine dorgbr( vect, m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  character :: vect
  integer :: info, k, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorgbr( vect, m, n, k, a, lda, tau, work, lwork, info )
end subroutine dorgbr

subroutine dorghr( n, ilo, ihi, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: ihi, ilo, info, lda, lwork, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorghr( n, ilo, ihi, a, lda, tau, work, lwork, info )
end subroutine dorghr

subroutine dorgl2( m, n, k, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, k, lda, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorgl2( m, n, k, a, lda, tau, work, info )
end subroutine dorgl2

subroutine dorglq( m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, k, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorglq( m, n, k, a, lda, tau, work, lwork, info )
end subroutine dorglq

subroutine dorgql( m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, k, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorgql( m, n, k, a, lda, tau, work, lwork, info )
end subroutine dorgql

subroutine dorgqr( m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, k, lda, lwork, m, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorgqr( m, n, k, a, lda, tau, work, lwork, info )
end subroutine dorgqr

subroutine dorgtr( uplo, n, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, lwork, n
  real*8 :: a( lda, * ), tau( * ), work( * )
  call lb_dorgtr( uplo, n, a, lda, tau, work, lwork, info )
end subroutine dorgtr

subroutine dorm2l( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dorm2l( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
end subroutine dorm2l

subroutine dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
end subroutine dorm2r

subroutine dormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans, vect
  integer :: info, k, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormbr

subroutine dormhr( side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: ihi, ilo, info, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormhr( side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormhr

subroutine dorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
end subroutine dorml2

subroutine dormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormlq

subroutine dormql( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormql( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormql

subroutine dormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans
  integer :: info, k, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormqr

subroutine dormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info )
  use link_blas
  implicit none
  character :: side, trans, uplo
  integer :: info, lda, ldc, lwork, m, n
  real*8 :: a( lda, * ), c( ldc, * ), tau( * ), work( * )
  call lb_dormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info )
end subroutine dormtr

subroutine dspev( jobz, uplo, n, ap, w, z, ldz, work, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer ::  info, ldz, n
  real*8 :: ap( * ), w( * ), work( * ), z( ldz, * )
  call lb_dspev( jobz, uplo, n, ap, w, z, ldz, work, info )
end subroutine dspev

subroutine dposv( uplo, n, nrhs, a, lda, b, ldb, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, ldb, n, nrhs
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dposv( uplo, n, nrhs, a, lda, b, ldb, info )
end subroutine dposv

subroutine dpotrf( uplo, n, a, lda, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, n
  real*8 :: a( lda, * )
  call lb_dpotrf( uplo, n, a, lda, info )
end subroutine dpotrf

recursive subroutine dpotrf2( uplo, n, a, lda, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, n
  real*8 :: a( lda, * )
  call lb_dpotrf2( uplo, n, a, lda, info )
end subroutine dpotrf2

subroutine dpotrs( uplo, n, nrhs, a, lda, b, ldb, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, ldb, n, nrhs
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dpotrs( uplo, n, nrhs, a, lda, b, ldb, info )
end subroutine dpotrs

subroutine dpptrf( uplo, n, ap, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, n
  real*8 :: ap( * )
  call lb_dpptrf( uplo, n, ap, info )
end subroutine dpptrf

subroutine drscl( n, sa, sx, incx )
  use link_blas
  implicit none
  integer :: incx, n
  real*8 :: sa
  real*8 :: sx( * )
  call lb_drscl( n, sa, sx, incx )
end subroutine drscl

subroutine dspgst( itype, uplo, n, ap, bp, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, itype, n
  real*8 :: ap( * ), bp( * )
  call lb_dspgst( itype, uplo, n, ap, bp, info )
end subroutine dspgst

subroutine dspgv( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer :: info, itype, ldz, n
  real*8 :: ap( * ), bp( * ), w( * ), work( * ), z( ldz, * )
  call lb_dspgv( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info )
end subroutine dspgv

subroutine dsptrd( uplo, n, ap, d, e, tau, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, n
  real*8 :: ap( * ), d( * ), e( * ), tau( * )
  call lb_dsptrd( uplo, n, ap, d, e, tau, info )
end subroutine dsptrd

subroutine dstebz( range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info )
  use link_blas
  implicit none
  character :: order, range
  integer :: il, info, iu, m, n, nsplit
  real*8 :: abstol, vl, vu
  integer :: iblock( * ), isplit( * ), iwork( * )
  real*8 :: d( * ), e( * ), w( * ), work( * )
  call lb_dstebz( range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info )
end subroutine dstebz

subroutine dstein( n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info )
  use link_blas
  implicit none
  integer :: info, ldz, m, n
  integer :: iblock( * ), ifail( * ), isplit( * ), iwork( * )
  real*8 :: d( * ), e( * ), w( * ), work( * ), z( ldz, * )
  call lb_dstein( n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info )
end subroutine dstein

subroutine dstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info )
  use link_blas
  implicit none
  character :: jobz, range
  logical :: tryrac
  integer :: il, info, iu, ldz, nzc, liwork, lwork, m, n
  real*8 :: vl, vu
  integer :: isuppz( * ), iwork( * )
  real*8 :: d( * ), e( * ), w( * ), work( * )
  real*8 :: z( ldz, * )
  call lb_dstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info )
end subroutine dstemr

subroutine dsteqr( compz, n, d, e, z, ldz, work, info )
  use link_blas
  implicit none
  character :: compz
  integer :: info, ldz, n
  real*8 :: d( * ), e( * ), work( * ), z( ldz, * )
  call lb_dsteqr( compz, n, d, e, z, ldz, work, info )
end subroutine dsteqr

subroutine dsterf( n, d, e, info )
  use link_blas
  implicit none
  integer :: info, n
  real*8 :: d( * ), e( * )
  call lb_dsterf( n, d, e, info )
end subroutine dsterf

subroutine dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer ::   info, lda, lwork, n
  real*8 ::    a( lda, * ), w( * ), work( * )
  call lb_dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
end subroutine dsyev

subroutine dstevr( jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info )
  use link_blas
  implicit none
  character :: jobz, range
  integer :: il, info, iu, ldz, liwork, lwork, m, n
  real*8 :: abstol, vl, vu
  integer :: isuppz( * ), iwork( * )
  real*8 :: d( * ), e( * ), w( * ), work( * ), z( ldz, * )
  call lb_dstevr( jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info )
end subroutine dstevr

subroutine dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info )
  use link_blas
  implicit none
  character :: jobz, range, uplo
  integer :: il, info, iu, lda, ldz, liwork, lwork, m, n
  real*8 :: abstol, vl, vu
  integer :: isuppz( * ), iwork( * )
  real*8 :: a( lda, * ), w( * ), work( * ), z( ldz, * )
  call lb_dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info )
end subroutine dsyevr

subroutine dsygs2( itype, uplo, n, a, lda, b, ldb, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, itype, lda, ldb, n
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dsygs2( itype, uplo, n, a, lda, b, ldb, info )
end subroutine dsygs2

subroutine dsygst( itype, uplo, n, a, lda, b, ldb, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, itype, lda, ldb, n
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dsygst( itype, uplo, n, a, lda, b, ldb, info )
end subroutine dsygst

subroutine dsygv( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer :: info, itype, lda, ldb, lwork, n
  real*8 :: a( lda, * ), b( ldb, * ), w( * ), work( * )
  call lb_dsygv( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info )
end subroutine dsygv

subroutine dsytd2( uplo, n, a, lda, d, e, tau, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, n
  real*8 :: a( lda, * ), d( * ), e( * ), tau( * )
  call lb_dsytd2( uplo, n, a, lda, d, e, tau, info )
end subroutine dsytd2

subroutine dsytrd( uplo, n, a, lda, d, e, tau, work, lwork, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, lwork, n
  real*8 :: a( lda, * ), d( * ), e( * ), tau( * ), work( * )
  call lb_dsytrd( uplo, n, a, lda, d, e, tau, work, lwork, info )
end subroutine dsytrd

subroutine dtrevc3( side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info )
  use link_blas
  implicit none
  character :: howmny, side
  integer :: info, ldt, ldvl, ldvr, lwork, m, mm, n
  logical :: select( * )
  real*8 :: t( ldt, * ), vl( ldvl, * ), vr( ldvr, * ), work( * )
  call lb_dtrevc3( side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info )
end subroutine dtrevc3

subroutine dtrexc( compq, n, t, ldt, q, ldq, ifst, ilst, work, info )
  use link_blas
  implicit none
  character :: compq
  integer :: ifst, ilst, info, ldq, ldt, n
  real*8 :: q( ldq, * ), t( ldt, * ), work( * )
  call lb_dtrexc( compq, n, t, ldt, q, ldq, ifst, ilst, work, info )
end subroutine dtrexc

subroutine dtrti2( uplo, diag, n, a, lda, info )
  use link_blas
  implicit none
  character :: diag, uplo
  integer :: info, lda, n
  real*8 :: a( lda, * )
  call lb_dtrti2( uplo, diag, n, a, lda, info )
end subroutine dtrti2

subroutine dtrtri( uplo, diag, n, a, lda, info )
  use link_blas
  implicit none
  character :: diag, uplo
  integer :: info, lda, n
  real*8 :: a( lda, * )
  call lb_dtrtri( uplo, diag, n, a, lda, info )
end subroutine dtrtri

subroutine dtrtrs( uplo, trans, diag, n, nrhs, a, lda, b, ldb, info )
  use link_blas
  implicit none
  character :: diag, trans, uplo
  integer :: info, lda, ldb, n, nrhs
  real*8 :: a( lda, * ), b( ldb, * )
  call lb_dtrtrs( uplo, trans, diag, n, nrhs, a, lda, b, ldb, info )
end subroutine dtrtrs

function ieeeck( ispec, zero, one )
  use link_blas
  implicit none
  integer :: ispec
  real :: one, zero
  integer :: ieeeck
  ieeeck=lb_ieeeck( ispec, zero, one )
end function ieeeck

function iladlc( m, n, a, lda )
  use link_blas
  implicit none
  integer :: m, n, lda
  real*8 :: a( lda, * )
  integer :: iladlc
  iladlc=lb_iladlc( m, n, a, lda )
end function iladlc

function iladlr( m, n, a, lda )
  use link_blas
  implicit none
  integer :: m, n, lda
  real*8 :: a( lda, * )
  integer :: iladlr
  iladlr=lb_iladlr( m, n, a, lda )
end function iladlr

function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
  use link_blas
  implicit none
  character*( * )    name, opts
  integer :: ispec, n1, n2, n3, n4
  integer :: ilaenv
  ilaenv=lb_ilaenv( ispec, name, opts, n1, n2, n3, n4 )
end function ilaenv

function ilazlc( m, n, a, lda )
  use link_blas
  implicit none
  integer :: m, n, lda
  complex*16 :: a( lda, * )
  integer :: ilazlc
  ilazlc=lb_ilazlc( m, n, a, lda )
end function ilazlc

function ilazlr( m, n, a, lda )
  use link_blas
  implicit none
  integer :: m, n, lda
  complex*16 :: a( lda, * )
  integer :: ilazlr
  ilazlr=lb_ilazlr( m, n, a, lda )
end function ilazlr

function iparmq( ispec, name, opts, n, ilo, ihi, lwork )
  use link_blas
  implicit none
  integer :: ihi, ilo, ispec, lwork, n
  character :: name*( * ), opts*( * )
  integer :: iparmq
  iparmq=lb_iparmq( ispec, name, opts, n, ilo, ihi, lwork )
end function iparmq

function iparam2stage( ispec, name, opts, ni, nbi, ibi, nxi )
  use link_blas
  implicit none
  character*( * ) :: name, opts
  integer :: ispec, ni, nbi, ibi, nxi
  integer :: iparam2stage
  iparam2stage=lb_iparam2stage( ispec, name, opts, ni, nbi, ibi, nxi )
end function iparam2stage

subroutine zheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer :: info, lda, lwork, n
  real*8 :: rwork( * ), w( * )
  complex*16 :: a( lda, * ), work( * )
  call lb_zheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info )
end subroutine zheev

subroutine zhetd2( uplo, n, a, lda, d, e, tau, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, n
  real*8 :: d( * ), e( * )
  complex*16 :: a( lda, * ), tau( * )
  call lb_zhetd2( uplo, n, a, lda, d, e, tau, info )
end subroutine zhetd2

subroutine zhetrd( uplo, n, a, lda, d, e, tau, work, lwork, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, lwork, n
  real*8 :: d( * ), e( * )
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zhetrd( uplo, n, a, lda, d, e, tau, work, lwork, info )
end subroutine zhetrd

subroutine zhpev( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
  use link_blas
  implicit none
  character :: jobz, uplo
  integer :: info, ldz, n
  real*8 :: rwork( * ), w( * )
  complex*16 :: ap( * ), work( * ), z( ldz, * )
  call lb_zhpev( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
end subroutine zhpev

subroutine zhptrd( uplo, n, ap, d, e, tau, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, n
  real*8 :: d( * ), e( * )
  complex*16 :: ap( * ), tau( * )
  call lb_zhptrd( uplo, n, ap, d, e, tau, info )
end subroutine zhptrd

subroutine zlacgv( n, x, incx )
  use link_blas
  implicit none
  integer :: incx, n
  complex*16 :: x( * )
  call lb_zlacgv( n, x, incx )
end subroutine zlacgv

function zladiv( x, y )
  use link_blas
  implicit none
  complex*16 :: x, y
  complex*16 :: zladiv
  zladiv=lb_zladiv( x, y )
end function zladiv

function zlanhe( norm, uplo, n, a, lda, work )
  use link_blas
  implicit none
  character :: norm, uplo
  integer :: lda, n
  real*8 :: work( * )
  complex*16 :: a( lda, * )
  real*8 :: zlanhe
  zlanhe=lb_zlanhe( norm, uplo, n, a, lda, work )
end function zlanhe

function zlanhp( norm, uplo, n, ap, work )
  use link_blas
  implicit none
  character :: norm, uplo
  integer :: n
  real*8 :: work( * )
  complex*16 :: ap( * )
  real*8 :: zlanhp
  zlanhp=lb_zlanhp( norm, uplo, n, ap, work )
end function zlanhp

subroutine zlarf( side, m, n, v, incv, tau, c, ldc, work )
  use link_blas
  implicit none
  character :: side
  integer :: incv, ldc, m, n
  complex*16 :: tau
  complex*16 :: c( ldc, * ), v( * ), work( * )
  call lb_zlarf( side, m, n, v, incv, tau, c, ldc, work )
end subroutine zlarf

subroutine zlarfb( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork )
  use link_blas
  implicit none
  character :: direct, side, storev, trans
  integer :: k, ldc, ldt, ldv, ldwork, m, n
  complex*16 :: c( ldc, * ), t( ldt, * ), v( ldv, * ), work( ldwork, * )
  call lb_zlarfb( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork )
end subroutine zlarfb

subroutine zlarfg( n, alpha, x, incx, tau )
  use link_blas
  implicit none
  integer :: incx, n
  complex*16 :: alpha, tau
  complex*16 :: x( * )
  call lb_zlarfg( n, alpha, x, incx, tau )
end subroutine zlarfg

subroutine zlarft( direct, storev, n, k, v, ldv, tau, t, ldt )
  use link_blas
  implicit none
  character :: direct, storev
  integer :: k, ldt, ldv, n
  complex*16 :: t( ldt, * ), tau( * ), v( ldv, * )
  call lb_zlarft( direct, storev, n, k, v, ldv, tau, t, ldt )
end subroutine zlarft

subroutine zlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
  use link_blas
  implicit none
  character :: type
  integer :: info, kl, ku, lda, m, n
  real*8 :: cfrom, cto
  complex*16 :: a( lda, * )
  call lb_zlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
end subroutine zlascl

subroutine zlaset( uplo, m, n, alpha, beta, a, lda )
  use link_blas
  implicit none
  character :: uplo
  integer :: lda, m, n
  complex*16 :: alpha, beta
  complex*16 :: a( lda, * )
  call lb_zlaset( uplo, m, n, alpha, beta, a, lda )
end subroutine zlaset

subroutine zlasr( side, pivot, direct, m, n, c, s, a, lda )
  use link_blas
  implicit none
  character :: direct, pivot, side
  integer :: lda, m, n
  real*8 :: c( * ), s( * )
  complex*16 :: a( lda, * )
  call lb_zlasr( side, pivot, direct, m, n, c, s, a, lda )
end subroutine zlasr

subroutine zlassq( n, x, incx, scale, sumsq )
  use link_blas
  implicit none
  integer :: incx, n
  real*8 :: scale, sumsq
  complex*16 :: x( * )
  call lb_zlassq( n, x, incx, scale, sumsq )
end subroutine zlassq

subroutine zlatrd( uplo, n, nb, a, lda, e, tau, w, ldw )
  use link_blas
  implicit none
  character :: uplo
  integer :: lda, ldw, n, nb
  real*8 :: e( * )
  complex*16 :: a( lda, * ), tau( * ), w( ldw, * )
  call lb_zlatrd( uplo, n, nb, a, lda, e, tau, w, ldw )
end subroutine zlatrd

subroutine zsteqr( compz, n, d, e, z, ldz, work, info )
  use link_blas
  implicit none
  character :: compz
  integer :: info, ldz, n
  real*8 :: d( * ), e( * ), work( * )
  complex*16 :: z( ldz, * )
  call lb_zsteqr( compz, n, d, e, z, ldz, work, info )
end subroutine zsteqr

subroutine zung2l( m, n, k, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, k, lda, m, n
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zung2l( m, n, k, a, lda, tau, work, info )
end subroutine zung2l

subroutine zung2r( m, n, k, a, lda, tau, work, info )
  use link_blas
  implicit none
  integer :: info, k, lda, m, n
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zung2r( m, n, k, a, lda, tau, work, info )
end subroutine zung2r

subroutine zungql( m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, k, lda, lwork, m, n
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zungql( m, n, k, a, lda, tau, work, lwork, info )
end subroutine zungql

subroutine zungqr( m, n, k, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  integer :: info, k, lda, lwork, m, n
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zungqr( m, n, k, a, lda, tau, work, lwork, info )
end subroutine zungqr

subroutine zungtr( uplo, n, a, lda, tau, work, lwork, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, lwork, n
  complex*16 :: a( lda, * ), tau( * ), work( * )
  call lb_zungtr( uplo, n, a, lda, tau, work, lwork, info )
end subroutine zungtr

subroutine zupgtr( uplo, n, ap, tau, q, ldq, work, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, ldq, n
  complex*16 :: ap( * ), q( ldq, * ), tau( * ), work( * )
  call lb_zupgtr( uplo, n, ap, tau, q, ldq, work, info )
end subroutine zupgtr

