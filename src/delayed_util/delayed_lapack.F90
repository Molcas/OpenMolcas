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

subroutine dbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork,info)
use link_blas, only: lb_dbdsdc
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compq, uplo
integer(kind=BLASInt) :: info, ldu, ldvt, n
integer(kind=BLASInt) :: iq(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), q(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork,info)
end subroutine dbdsdc

subroutine dbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)
use link_blas, only: lb_dbdsqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, ldc, ldu, ldvt, n, ncc, ncvt, nru
real(kind=BLASR8) :: c(ldc,*), d(*), e(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)
end subroutine dbdsqr

subroutine dgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
use link_blas, only: lb_dgebak
use Definitions, only: BLASInt, BLASR8
implicit none
character :: job, side
integer(kind=BLASInt) :: ihi, ilo, info, ldv, m, n
real(kind=BLASR8) :: scale(*), v(ldv,*)
call lb_dgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
end subroutine dgebak

subroutine dgebal(job,n,a,lda,ilo,ihi,scale,info)
use link_blas, only: lb_dgebal
use Definitions, only: BLASInt, BLASR8
implicit none
character :: job
integer(kind=BLASInt) :: ihi, ilo, info, lda, n
real(kind=BLASR8) :: a(lda,*), scale(*)
call lb_dgebal(job,n,a,lda,ilo,ihi,scale,info)
end subroutine dgebal

subroutine dgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
use link_blas, only: lb_dgebd2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
real(kind=BLASR8) :: a(lda,*), d(*), e(*), taup(*), tauq(*), work(*)
call lb_dgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
end subroutine dgebd2

subroutine dgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
use link_blas, only: lb_dgebrd
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), d(*), e(*), taup(*), tauq(*), work(*)
call lb_dgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
end subroutine dgebrd

subroutine dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)
use link_blas, only: lb_dgecon
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: anorm, rcond
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: a(lda,*), work(*)
call lb_dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)
end subroutine dgecon

subroutine dgees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork,bwork,info)
use link_blas, only: lb_dgees
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobvs, sort
integer(kind=BLASInt) :: info, lda, ldvs, lwork, n, sdim
logical :: bwork(*)
real(kind=BLASR8) :: a(lda,*), vs(ldvs,*), wi(*), work(*), wr(*)
logical, external :: select
call lb_dgees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork,bwork,info)
end subroutine dgees

subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
use link_blas, only: lb_dgeev
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobvl, jobvr
integer(kind=BLASInt) :: info, lda, ldvl, ldvr, lwork, n
real(kind=BLASR8) :: a(lda,*), vl(ldvl,*), vr(ldvr,*), wi(*), work(*), wr(*)
call lb_dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
end subroutine dgeev

subroutine dgehd2(n,ilo,ihi,a,lda,tau,work,info)
use link_blas, only: lb_dgehd2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ilo, info, lda, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgehd2(n,ilo,ihi,a,lda,tau,work,info)
end subroutine dgehd2

subroutine dgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dgehrd
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ilo, info, lda, lwork, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
end subroutine dgehrd

subroutine dgelq2(m,n,a,lda,tau,work,info)
use link_blas, only: lb_dgelq2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgelq2(m,n,a,lda,tau,work,info)
end subroutine dgelq2

subroutine dgelqf(m,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dgelqf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgelqf(m,n,a,lda,tau,work,lwork,info)
end subroutine dgelqf

subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
use link_blas, only: lb_dgels
use Definitions, only: BLASInt, BLASR8
implicit none
character :: trans
integer(kind=BLASInt) :: info, lda, ldb, lwork, m, n, nrhs
real(kind=BLASR8) :: a(lda,*), b(ldb,*), work(*)
call lb_dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
end subroutine dgels

subroutine dgeqr2(m,n,a,lda,tau,work,info)
use link_blas, only: lb_dgeqr2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgeqr2(m,n,a,lda,tau,work,info)
end subroutine dgeqr2

subroutine dgeqrf(m,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dgeqrf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dgeqrf(m,n,a,lda,tau,work,lwork,info)
end subroutine dgeqrf

subroutine dgesdd(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)
use link_blas, only: lb_dgesdd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz
integer(kind=BLASInt) :: info, lda, ldu, ldvt, lwork, m, n
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dgesdd(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)
end subroutine dgesdd

subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
use link_blas, only: lb_dgesv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
end subroutine dgesv

subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
use link_blas, only: lb_dgesvd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobu, jobvt
integer(kind=BLASInt) :: info, lda, ldu, ldvt, lwork, m, n
real(kind=BLASR8) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
end subroutine dgesvd

subroutine dgetrf(m,n,a,lda,ipiv,info)
use link_blas, only: lb_dgetrf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*)
call lb_dgetrf(m,n,a,lda,ipiv,info)
end subroutine dgetrf

recursive subroutine dgetrf2(m,n,a,lda,ipiv,info)
use link_blas, only: lb_dgetrf2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*)
call lb_dgetrf2(m,n,a,lda,ipiv,info)
end subroutine dgetrf2

subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
use link_blas, only: lb_dgetri
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, n
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*), work(*)
call lb_dgetri(n,a,lda,ipiv,work,lwork,info)
end subroutine dgetri

subroutine dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
use link_blas, only: lb_dgetrs
use Definitions, only: BLASInt, BLASR8
implicit none
character :: trans
integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
end subroutine dgetrs

subroutine dhseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
use link_blas, only: lb_dhseqr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ilo, info, ldh, ldz, lwork, n
character :: compz, job
real(kind=BLASR8) :: h(ldh,*), wi(*), work(*), wr(*), z(ldz,*)
call lb_dhseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
end subroutine dhseqr

function disnan(din)
use link_blas, only: lb_disnan
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: din
logical :: disnan
disnan = lb_disnan(din)
end function disnan

subroutine dlabad(small,large)
use link_blas, only: lb_dlabad
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: large, small
call lb_dlabad(small,large)
end subroutine dlabad

subroutine dlabrd(m,n,nb,a,lda,d,e,tauq,taup,x,ldx,y,ldy)
use link_blas, only: lb_dlabrd
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: lda, ldx, ldy, m, n, nb
real(kind=BLASR8) :: a(lda,*), d(*), e(*), taup(*), tauq(*), x(ldx,*), y(ldy,*)
call lb_dlabrd(m,n,nb,a,lda,d,e,tauq,taup,x,ldx,y,ldy)
end subroutine dlabrd

subroutine dlacn2(n,v,x,isgn,est,kase,isave)
use link_blas, only: lb_dlacn2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: kase, n
real(kind=BLASR8) :: est
integer(kind=BLASInt) :: isgn(*), isave(3)
real(kind=BLASR8) :: v(*), x(*)
call lb_dlacn2(n,v,x,isgn,est,kase,isave)
end subroutine dlacn2

subroutine dlacpy(uplo,m,n,a,lda,b,ldb)
use link_blas, only: lb_dlacpy
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, ldb, m, n
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dlacpy(uplo,m,n,a,lda,b,ldb)
end subroutine dlacpy

subroutine dladiv(a,b,c,d,p,q)
use link_blas, only: lb_dladiv
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, d, p, q
call lb_dladiv(a,b,c,d,p,q)
end subroutine dladiv

subroutine dladiv1(a,b,c,d,p,q)
use link_blas, only: lb_dladiv1
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, d, p, q
call lb_dladiv1(a,b,c,d,p,q)
end subroutine dladiv1

function dladiv2(a,b,c,d,r,t)
use link_blas, only: lb_dladiv2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, d, r, t
real(kind=BLASR8) :: dladiv2
dladiv2 = lb_dladiv2(a,b,c,d,r,t)
end function dladiv2

subroutine dlae2(a,b,c,rt1,rt2)
use link_blas, only: lb_dlae2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, rt1, rt2
call lb_dlae2(a,b,c,rt1,rt2)
end subroutine dlae2

subroutine dlaebz(ijob,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,e,e2,nval,ab,c,mout,nab,work,iwork,info)
use link_blas, only: lb_dlaebz
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ijob, info, minp, mmax, mout, n, nbmin, nitmax
real(kind=BLASR8) :: abstol, pivmin, reltol
integer(kind=BLASInt) :: iwork(*), nab(mmax,*), nval(*)
real(kind=BLASR8) :: ab(mmax,*), c(*), d(*), e(*), e2(*), work(*)
call lb_dlaebz(ijob,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,e,e2,nval,ab,c,mout,nab,work,iwork,info)
end subroutine dlaebz

subroutine dlaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info)
use link_blas, only: lb_dlaed0
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: icompq, info, ldq, ldqs, n, qsiz
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e(*), q(ldq,*), qstore(ldqs,*), work(*)
call lb_dlaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info)
end subroutine dlaed0

subroutine dlaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
use link_blas, only: lb_dlaed1
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: cutpnt, info, ldq, n
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: indxq(*), iwork(*)
real(kind=BLASR8) :: d(*), q(ldq,*), work(*)
call lb_dlaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
end subroutine dlaed1

subroutine dlaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc,indxp,coltyp,info)
use link_blas, only: lb_dlaed2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, ldq, n, n1
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: coltyp(*), indx(*), indxc(*), indxp(*), indxq(*)
real(kind=BLASR8) :: d(*), dlamda(*), q(ldq,*), q2(*), w(*), z(*)
call lb_dlaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc,indxp,coltyp,info)
end subroutine dlaed2

subroutine dlaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
use link_blas, only: lb_dlaed3
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, ldq, n, n1
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: ctot(*), indx(*)
real(kind=BLASR8) :: d(*), dlamda(*), q(ldq,*), q2(*), s(*), w(*)
call lb_dlaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
end subroutine dlaed3

subroutine dlaed4(n,i,d,z,delta,rho,dlam,info)
use link_blas, only: lb_dlaed4
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i, info, n
real(kind=BLASR8) :: dlam, rho
real(kind=BLASR8) :: d(*), delta(*), z(*)
call lb_dlaed4(n,i,d,z,delta,rho,dlam,info)
end subroutine dlaed4

subroutine dlaed5(i,d,z,delta,rho,dlam)
use link_blas, only: lb_dlaed5
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i
real(kind=BLASR8) :: dlam, rho
real(kind=BLASR8) :: d(2), delta(2), z(2)
call lb_dlaed5(i,d,z,delta,rho,dlam)
end subroutine dlaed5

subroutine dlaed6(kniter,orgati,rho,d,z,finit,tau,info)
use link_blas, only: lb_dlaed6
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: orgati
integer(kind=BLASInt) :: info, kniter
real(kind=BLASR8) :: finit, rho, tau
real(kind=BLASR8) :: d(3), z(3)
call lb_dlaed6(kniter,orgati,rho,d,z,finit,tau,info)
end subroutine dlaed6

subroutine dlaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho,cutpnt,qstore,qptr,prmptr,perm,givptr, &
                  givcol,givnum,work,iwork,info)
use link_blas, only: lb_dlaed7
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: curlvl, curpbm, cutpnt, icompq, info, ldq, n, qsiz, tlvls
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: givcol(2,*), givptr(*), indxq(*), iwork(*), perm(*), prmptr(*), qptr(*)
real(kind=BLASR8) :: d(*), givnum(2,*), q(ldq,*), qstore(*), work(*)
call lb_dlaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho,cutpnt,qstore,qptr,prmptr,perm,givptr, &
               givcol,givnum,work,iwork,info)
end subroutine dlaed7

subroutine dlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda,q2,ldq2,w,perm,givptr,givcol,givnum,indxp, &
                  indx,info)
use link_blas, only: lb_dlaed8
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: cutpnt, givptr, icompq, info, k, ldq, ldq2, n, qsiz
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: givcol(2,*), indx(*), indxp(*), indxq(*), perm(*)
real(kind=BLASR8) :: d(*), dlamda(*), givnum(2,*), q(ldq,*), q2(ldq2,*), w(*), z(*)
call lb_dlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda,q2,ldq2,w,perm,givptr,givcol,givnum,indxp, &
               indx,info)
end subroutine dlaed8

subroutine dlaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)
use link_blas, only: lb_dlaed9
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, kstart, kstop, ldq, lds, n
real(kind=BLASR8) :: rho
real(kind=BLASR8) :: d(*), dlamda(*), q(ldq,*), s(lds,*), w(*)
call lb_dlaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)
end subroutine dlaed9

subroutine dlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum,q,qptr,z,ztemp,info)
use link_blas, only: lb_dlaeda
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: curlvl, curpbm, info, n, tlvls
integer(kind=BLASInt) :: givcol(2,*), givptr(*), perm(*), prmptr(*), qptr(*)
real(kind=BLASR8) :: givnum(2,*), q(*), z(*), ztemp(*)
call lb_dlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum,q,qptr,z,ztemp,info)
end subroutine dlaeda

subroutine dlaev2(a,b,c,rt1,rt2,cs1,sn1)
use link_blas, only: lb_dlaev2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, cs1, rt1, rt2, sn1
call lb_dlaev2(a,b,c,rt1,rt2,cs1,sn1)
end subroutine dlaev2

subroutine dlaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
use link_blas, only: lb_dlaexc
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: wantq
integer(kind=BLASInt) :: info, j1, ldq, ldt, n, n1, n2
real(kind=BLASR8) :: q(ldq,*), t(ldt,*), work(*)
call lb_dlaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
end subroutine dlaexc

subroutine dlagtf(n,a,lambda,b,c,tol,d,in,info)
use link_blas, only: lb_dlagtf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: lambda, tol
integer(kind=BLASInt) :: in(*)
real(kind=BLASR8) :: a(*), b(*), c(*), d(*)
call lb_dlagtf(n,a,lambda,b,c,tol,d,in,info)
end subroutine dlagtf

subroutine dlagts(job,n,a,b,c,d,in,y,tol,info)
use link_blas, only: lb_dlagts
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, job, n
real(kind=BLASR8) :: tol
integer(kind=BLASInt) :: in(*)
real(kind=BLASR8) :: a(*), b(*), c(*), d(*), y(*)
call lb_dlagts(job,n,a,b,c,d,in,y,tol,info)
end subroutine dlagts

subroutine dlahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,info)
use link_blas, only: lb_dlahqr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ihiz, ilo, iloz, info, ldh, ldz, n
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), wi(*), wr(*), z(ldz,*)
call lb_dlahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,info)
end subroutine dlahqr

subroutine dlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
use link_blas, only: lb_dlahr2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: k, lda, ldt, ldy, n, nb
real(kind=BLASR8) :: a(lda,*), t(ldt,nb), tau(nb), y(ldy,nb)
call lb_dlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
end subroutine dlahr2

function dlaisnan(din1,din2)
use link_blas, only: lb_dlaisnan
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: din1, din2
logical :: dlaisnan
dlaisnan = lb_dlaisnan(din1,din2)
end function dlaisnan

subroutine dlaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x,ldx,scale,xnorm,info)
use link_blas, only: lb_dlaln2
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: ltrans
integer(kind=BLASInt) :: info, lda, ldb, ldx, na, nw
real(kind=BLASR8) :: ca, d1, d2, scale, smin, wi, wr, xnorm
real(kind=BLASR8) :: a(lda,*), b(ldb,*), x(ldx,*)
call lb_dlaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x,ldx,scale,xnorm,info)
end subroutine dlaln2

function dlamc3(a,b)
use link_blas, only: lb_dlamc3
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: A, B
real(kind=BLASR8) :: dlamc3
dlamc3 = lb_dlamc3(a,b)
end function dlamc3

function dlamch(cmach)
use link_blas, only: lb_dlamch
use Definitions, only: BLASR8
implicit none
character :: cmach
real(kind=BLASR8) :: dlamch
dlamch = lb_dlamch(cmach)
end function dlamch

subroutine dlamrg(n1,n2,a,dtrd1,dtrd2,index)
use link_blas, only: lb_dlamrg
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: dtrd1, dtrd2, n1, n2
integer(kind=BLASInt) :: index(*)
real(kind=BLASR8) :: a(*)
call lb_dlamrg(n1,n2,a,dtrd1,dtrd2,index)
end subroutine dlamrg

function dlaneg(n,d,lld,sigma,pivmin,r)
use link_blas, only: lb_dlaneg
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: n, r
real(kind=BLASR8) :: pivmin, sigma
real(kind=BLASR8) :: d(*), lld(*)
integer(kind=BLASInt) :: dlaneg
dlaneg = lb_dlaneg(n,d,lld,sigma,pivmin,r)
end function dlaneg

function dlange(norm,m,n,a,lda,work)
use link_blas, only: lb_dlange
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm
integer(kind=BLASInt) :: lda, m, n
real(kind=BLASR8) :: a(lda,*), work(*)
real(kind=BLASR8) :: dlange
dlange = lb_dlange(norm,m,n,a,lda,work)
end function dlange

function dlansp(norm,uplo,n,ap,work)
use link_blas, only: lb_dlansp
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm, uplo
integer(kind=BLASInt) :: n
real(kind=BLASR8) :: ap(*), work(*)
real(kind=BLASR8) :: dlansp
dlansp = lb_dlansp(norm,uplo,n,ap,work)
end function dlansp

function dlanst(norm,n,d,e)
use link_blas, only: lb_dlanst
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm
integer(kind=BLASInt) :: n
real(kind=BLASR8) :: d(*), e(*)
real(kind=BLASR8) :: dlanst
dlanst = lb_dlanst(norm,n,d,e)
end function dlanst

function dlansy(norm,uplo,n,a,lda,work)
use link_blas, only: lb_dlansy
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm, uplo
integer(kind=BLASInt) :: lda, n
real(kind=BLASR8) :: a(lda,*), work(*)
real(kind=BLASR8) :: dlansy
dlansy = lb_dlansy(norm,uplo,n,a,lda,work)
end function dlansy

subroutine dlanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
use link_blas, only: lb_dlanv2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: a, b, c, cs, d, rt1i, rt1r, rt2i, rt2r, sn
call lb_dlanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
end subroutine dlanv2

function dlapy2(x,y)
use link_blas, only: lb_dlapy2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: x, y
real(kind=BLASR8) :: dlapy2
dlapy2 = lb_dlapy2(x,y)
end function dlapy2

function dlapy3(x,y,z)
use link_blas, only: lb_dlapy3
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: x, y, z
real(kind=BLASR8) :: dlapy3
dlapy3 = lb_dlapy3(x,y,z)
end function dlapy3

subroutine dlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,work,lwork,info)
use link_blas, only: lb_dlaqr0
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ihiz, ilo, iloz, info, ldh, ldz, lwork, n
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), wi(*), work(*), wr(*), z(ldz,*)
call lb_dlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,work,lwork,info)
end subroutine dlaqr0

subroutine dlaqr1(n,h,ldh,sr1,si1,sr2,si2,v)
use link_blas, only: lb_dlaqr1
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: si1, si2, sr1, sr2
integer(kind=BLASInt) :: ldh, n
real(kind=BLASR8) :: h(ldh,*), v(*)
call lb_dlaqr1(n,h,ldh,sr1,si1,sr2,si2,v)
end subroutine dlaqr1

subroutine dlaqr2(wantt,wantz,n,ktop,kbot,nw,h,ldh,iloz,ihiz,z,ldz,ns,nd,sr,si,v,ldv,nh,t,ldt,nv,wv,ldwv, &
                  work,lwork)
use link_blas, only: lb_dlaqr2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihiz, iloz, kbot, ktop, ldh, ldt, ldv, ldwv, ldz, lwork, n, nd, nh, ns, nv, nw
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), si(*), sr(*), t(ldt,*), v(ldv,*), work(*), wv(ldwv,*), z(ldz,*)
call lb_dlaqr2(wantt,wantz,n,ktop,kbot,nw,h,ldh,iloz,ihiz,z,ldz,ns,nd,sr,si,v,ldv,nh,t,ldt,nv,wv,ldwv, &
               work,lwork)
end subroutine dlaqr2

subroutine dlaqr3(wantt,wantz,n,ktop,kbot,nw,h,ldh,iloz,ihiz,z,ldz,ns,nd,sr,si,v,ldv,nh,t,ldt,nv,wv,ldwv, &
                  work,lwork)
use link_blas, only: lb_dlaqr3
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihiz, iloz, kbot, ktop, ldh, ldt, ldv, ldwv, ldz, lwork, n, nd, nh, ns, nv, nw
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), si(*), sr(*), t(ldt,*), v(ldv,*), work(*), wv(ldwv,*), z(ldz,*)
call lb_dlaqr3(wantt,wantz,n,ktop,kbot,nw,h,ldh,iloz,ihiz,z,ldz,ns,nd,sr,si,v,ldv,nh,t,ldt,nv,wv,ldwv, &
               work,lwork)
end subroutine dlaqr3

subroutine dlaqr4(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,work,lwork,info)
use link_blas, only: lb_dlaqr4
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ihiz, ilo, iloz, info, ldh, ldz, lwork, n
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), wi(*), work(*), wr(*), z(ldz,*)
call lb_dlaqr4(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,work,lwork,info)
end subroutine dlaqr4

subroutine dlaqr5(wantt,wantz,kacc22,n,ktop,kbot,nshfts,sr,si,h,ldh,iloz,ihiz,z,ldz,v,ldv,u,ldu,nv,wv,ldwv, &
                  nh,wh,ldwh)
use link_blas, only: lb_dlaqr5
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihiz, iloz, kacc22, kbot, ktop, ldh, ldu, ldv, ldwh, ldwv, ldz, n, nh, nshfts, nv
logical :: wantt, wantz
real(kind=BLASR8) :: h(ldh,*), si(*), sr(*), u(ldu,*), v(ldv,*), wh(ldwh,*), wv(ldwv,*), z(ldz,*)
call lb_dlaqr5(wantt,wantz,kacc22,n,ktop,kbot,nshfts,sr,si,h,ldh,iloz,ihiz,z,ldz,v,ldv,u,ldu,nv,wv,ldwv, &
               nh,wh,ldwh)
end subroutine dlaqr5

subroutine dlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc,negcnt,ztz,mingma,r,isuppz,nrminv,resid, &
                  rqcorr,work)
use link_blas, only: lb_dlar1v
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: wantnc
integer(kind=BLASInt) :: b1, bn, n, negcnt, r
real(kind=BLASR8) :: gaptol, lambda, mingma, nrminv, pivmin, resid, rqcorr, ztz
integer(kind=BLASInt) :: isuppz(*)
real(kind=BLASR8) :: d(*), l(*), ld(*), lld(*), work(*)
real(kind=BLASR8) :: z(*)
call lb_dlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc,negcnt,ztz,mingma,r,isuppz,nrminv,resid, &
               rqcorr,work)
end subroutine dlar1v

subroutine dlarf(side,m,n,v,incv,tau,c,ldc,work)
use link_blas, only: lb_dlarf
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side
integer(kind=BLASInt) :: incv, ldc, m, n
real(kind=BLASR8) :: tau
real(kind=BLASR8) :: c(ldc,*), v(*), work(*)
call lb_dlarf(side,m,n,v,incv,tau,c,ldc,work)
end subroutine dlarf

subroutine dlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc,work,ldwork)
use link_blas, only: lb_dlarfb
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, side, storev, trans
integer(kind=BLASInt) :: k, ldc, ldt, ldv, ldwork, m, n
real(kind=BLASR8) :: c(ldc,*), t(ldt,*), v(ldv,*), work(ldwork,*)
call lb_dlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc,work,ldwork)
end subroutine dlarfb

subroutine dlarfg(n,alpha,x,incx,tau)
use link_blas, only: lb_dlarfg
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: alpha, tau
real(kind=BLASR8) :: x(*)
call lb_dlarfg(n,alpha,x,incx,tau)
end subroutine dlarfg

subroutine dlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
use link_blas, only: lb_dlarft
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, storev
integer(kind=BLASInt) :: k, ldt, ldv, n
real(kind=BLASR8) :: t(ldt,*), tau(*), v(ldv,*)
call lb_dlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
end subroutine dlarft

subroutine dlarfx(side,m,n,v,tau,c,ldc,work)
use link_blas, only: lb_dlarfx
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side
integer(kind=BLASInt) :: ldc, m, n
real(kind=BLASR8) :: tau
real(kind=BLASR8) :: c(ldc,*), v(*), work(*)
call lb_dlarfx(side,m,n,v,tau,c,ldc,work)
end subroutine dlarfx

subroutine dlarnv(idist,iseed,n,x)
use link_blas, only: lb_dlarnv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: idist, n
integer(kind=BLASInt) :: iseed(4)
real(kind=BLASR8) :: x(*)
call lb_dlarnv(idist,iseed,n,x)
end subroutine dlarnv

subroutine dlarra(n,d,e,e2,spltol,tnrm,nsplit,isplit,info)
use link_blas, only: lb_dlarra
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, n, nsplit
real(kind=BLASR8) :: spltol, tnrm
integer(kind=BLASInt) :: isplit(*)
real(kind=BLASR8) :: d(*), e(*), e2(*)
call lb_dlarra(n,d,e,e2,spltol,tnrm,nsplit,isplit,info)
end subroutine dlarra

subroutine dlarrb(n,d,lld,ifirst,ilast,rtol1,rtol2,offset,w,wgap,werr,work,iwork,pivmin,spdiam,twist,info)
use link_blas, only: lb_dlarrb
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ifirst, ilast, info, n, offset, twist
real(kind=BLASR8) :: pivmin, rtol1, rtol2, spdiam
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), lld(*), w(*), werr(*), wgap(*), work(*)
call lb_dlarrb(n,d,lld,ifirst,ilast,rtol1,rtol2,offset,w,wgap,werr,work,iwork,pivmin,spdiam,twist,info)
end subroutine dlarrb

subroutine dlarrc(jobt,n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info)
use link_blas, only: lb_dlarrc
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobt
integer(kind=BLASInt) :: eigcnt, info, lcnt, n, rcnt
real(kind=BLASR8) :: pivmin, vl, vu
real(kind=BLASR8) :: d(*), e(*)
call lb_dlarrc(jobt,n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info)
end subroutine dlarrc

subroutine dlarrd(range,order,n,vl,vu,il,iu,gers,reltol,d,e,e2,pivmin,nsplit,isplit,m,w,werr,wl,wu,iblock, &
                  indexw,work,iwork,info)
use link_blas, only: lb_dlarrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: order, range
integer(kind=BLASInt) :: il, info, iu, m, n, nsplit
real(kind=BLASR8) :: pivmin, reltol, vl, vu, wl, wu
integer(kind=BLASInt) :: iblock(*), indexw(*), isplit(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), e2(*), gers(*), w(*), werr(*), work(*)
call lb_dlarrd(range,order,n,vl,vu,il,iu,gers,reltol,d,e,e2,pivmin,nsplit,isplit,m,w,werr,wl,wu,iblock, &
               indexw,work,iwork,info)
end subroutine dlarrd

subroutine dlarre(range,n,vl,vu,il,iu,d,e,e2,rtol1,rtol2,spltol,nsplit,isplit,m,w,werr,wgap,iblock,indexw, &
                  gers,pivmin,work,iwork,info)
use link_blas, only: lb_dlarre
use Definitions, only: BLASInt, BLASR8
implicit none
character :: range
integer(kind=BLASInt) :: il, info, iu, m, n, nsplit
real(kind=BLASR8) :: pivmin, rtol1, rtol2, spltol, vl, vu
integer(kind=BLASInt) :: iblock(*), isplit(*), iwork(*), indexw(*)
real(kind=BLASR8) :: d(*), e(*), e2(*), gers(*), w(*), werr(*), wgap(*), work(*)
call lb_dlarre(range,n,vl,vu,il,iu,d,e,e2,rtol1,rtol2,spltol,nsplit,isplit,m,w,werr,wgap,iblock,indexw, &
               gers,pivmin,work,iwork,info)
end subroutine dlarre

subroutine dlarrf(n,d,l,ld,clstrt,clend,w,wgap,werr,spdiam,clgapl,clgapr,pivmin,sigma,dplus,lplus,work,info)
use link_blas, only: lb_dlarrf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: clstrt, clend, info, n
real(kind=BLASR8) :: clgapl, clgapr, pivmin, sigma, spdiam
real(kind=BLASR8) :: d(*), dplus(*), l(*), ld(*), lplus(*), w(*), wgap(*), werr(*), work(*)
call lb_dlarrf(n,d,l,ld,clstrt,clend,w,wgap,werr,spdiam,clgapl,clgapr,pivmin,sigma,dplus,lplus,work,info)
end subroutine dlarrf

subroutine dlarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork,pivmin,spdiam,info)
use link_blas, only: lb_dlarrj
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ifirst, ilast, info, n, offset
real(kind=BLASR8) :: pivmin, rtol, spdiam
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e2(*), w(*), werr(*), work(*)
call lb_dlarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork,pivmin,spdiam,info)
end subroutine dlarrj

subroutine dlarrk(n,iw,gl,gu,d,e2,pivmin,reltol,w,werr,info)
use link_blas, only: lb_dlarrk
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, iw, n
real(kind=BLASR8) :: pivmin, reltol, gl, gu, w, werr
real(kind=BLASR8) :: d(*), e2(*)
call lb_dlarrk(n,iw,gl,gu,d,e2,pivmin,reltol,w,werr,info)
end subroutine dlarrk

subroutine dlarrr(n,d,e,info)
use link_blas, only: lb_dlarrr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: n, info
real(kind=BLASR8) :: d(*), e(*)
call lb_dlarrr(n,d,e,info)
end subroutine dlarrr

subroutine dlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1,rtol2,w,werr,wgap,iblock,indexw,gers,z, &
                  ldz,isuppz,work,iwork,info)
use link_blas, only: lb_dlarrv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: dol, dou, info, ldz, m, n
real(kind=BLASR8) :: minrgp, pivmin, rtol1, rtol2, vl, vu
integer(kind=BLASInt) :: iblock(*), indexw(*), isplit(*), isuppz(*), iwork(*)
real(kind=BLASR8) :: d(*), gers(*), l(*), w(*), werr(*), wgap(*), work(*)
real(kind=BLASR8) :: z(ldz,*)
call lb_dlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1,rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz, &
               isuppz,work,iwork,info)
end subroutine dlarrv

subroutine dlartg(f,g,cs,sn,r)
use link_blas, only: lb_dlartg
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: cs, f, g, r, sn
call lb_dlartg(f,g,cs,sn,r)
end subroutine dlartg

subroutine dlaruv(iseed,n,x)
use link_blas, only: lb_dlaruv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: n
integer(kind=BLASInt) :: iseed(4)
real(kind=BLASR8) :: x(n)
call lb_dlaruv(iseed,n,x)
end subroutine dlaruv

subroutine dlas2(f,g,h,ssmin,ssmax)
use link_blas, only: lb_dlas2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: f, g, h, ssmax, ssmin
call lb_dlas2(f,g,h,ssmin,ssmax)
end subroutine dlas2

subroutine dlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
use link_blas, only: lb_dlascl
use Definitions, only: BLASInt, BLASR8
implicit none
character :: type
integer(kind=BLASInt) :: info, kl, ku, lda, m, n
real(kind=BLASR8) :: cfrom, cto
real(kind=BLASR8) :: a(lda,*)
call lb_dlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
end subroutine dlascl

subroutine dlasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)
use link_blas, only: lb_dlasd0
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, ldu, ldvt, n, smlsiz, sqre
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dlasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)
end subroutine dlasd0

subroutine dlasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork,work,info)
use link_blas, only: lb_dlasd1
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, ldu, ldvt, nl, nr, sqre
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: idxq(*), iwork(*)
real(kind=BLASR8) :: d(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dlasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork,work,info)
end subroutine dlasd1

subroutine dlasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma,u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq, &
                  coltyp,info)
use link_blas, only: lb_dlasd2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, ldu, ldu2, ldvt, ldvt2, nl, nr, sqre
real(kind=BLASR8) :: alpha, beta
integer(kind=BLASInt) :: coltyp(*), idx(*), idxc(*), idxp(*), idxq(*)
real(kind=BLASR8) :: d(*), dsigma(*), u(ldu,*), u2(ldu2,*), vt(ldvt,*), vt2(ldvt2,*), z(*)
call lb_dlasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma,u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq, &
               coltyp,info)
end subroutine dlasd2

subroutine dlasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt,vt2,ldvt2,idxc,ctot,z,info)
use link_blas, only: lb_dlasd3
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, ldq, ldu, ldu2, ldvt, ldvt2, nl, nr, sqre
integer(kind=BLASInt) :: ctot(*), idxc(*)
real(kind=BLASR8) :: d(*), dsigma(*), q(ldq,*), u(ldu,*), u2(ldu2,*), vt(ldvt,*), vt2(ldvt2,*), z(*)
call lb_dlasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt,vt2,ldvt2,idxc,ctot,z,info)
end subroutine dlasd3

subroutine dlasd4(n,i,d,z,delta,rho,sigma,work,info)
use link_blas, only: lb_dlasd4
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i, info, n
real(kind=BLASR8) :: rho, sigma
real(kind=BLASR8) :: d(*), delta(*), work(*), z(*)
call lb_dlasd4(n,i,d,z,delta,rho,sigma,work,info)
end subroutine dlasd4

subroutine dlasd5(i,d,z,delta,rho,dsigma,work)
use link_blas, only: lb_dlasd5
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i
real(kind=BLASR8) :: dsigma, rho
real(kind=BLASR8) :: d(2), delta(2), work(2), z(2)
call lb_dlasd5(i,d,z,delta,rho,dsigma,work)
end subroutine dlasd5

subroutine dlasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,poles,difl, &
                  difr,z,k,c,s,work,iwork,info)
use link_blas, only: lb_dlasd6
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: givptr, icompq, info, k, ldgcol, ldgnum, nl, nr, sqre
real(kind=BLASR8) :: alpha, beta, c, s
integer(kind=BLASInt) :: givcol(ldgcol,*), idxq(*), iwork(*), perm(*)
real(kind=BLASR8) :: d(*), difl(*), difr(*), givnum(ldgnum,*), poles(ldgnum,*), vf(*), vl(*), work(*), z(*)
call lb_dlasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,poles,difl, &
               difr,z,k,c,s,work,iwork,info)
end subroutine dlasd6

subroutine dlasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha,beta,dsigma,idx,idxp,idxq,perm,givptr, &
                  givcol,ldgcol,givnum,ldgnum,c,s,info)
use link_blas, only: lb_dlasd7
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: givptr, icompq, info, k, ldgcol, ldgnum, nl, nr, sqre
real(kind=BLASR8) :: alpha, beta, c, s
integer(kind=BLASInt) :: givcol(ldgcol,*), idx(*), idxp(*), idxq(*), perm(*)
real(kind=BLASR8) :: d(*), dsigma(*), givnum(ldgnum,*), vf(*), vfw(*), vl(*), vlw(*), z(*), zw(*)
call lb_dlasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha,beta,dsigma,idx,idxp,idxq,perm,givptr, &
               givcol,ldgcol,givnum,ldgnum,c,s,info)
end subroutine dlasd7

subroutine dlasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work,info)
use link_blas, only: lb_dlasd8
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: icompq, info, k, lddifr
real(kind=BLASR8) :: d(*), difl(*), difr(lddifr,*), dsigma(*), vf(*), vl(*), work(*), z(*)
call lb_dlasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work,info)
end subroutine dlasd8

subroutine dlasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z,poles,givptr,givcol,ldgcol,perm,givnum,c,s, &
                  work,iwork,info)
use link_blas, only: lb_dlasda
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: icompq, info, ldgcol, ldu, n, smlsiz, sqre
integer(kind=BLASInt) :: givcol(ldgcol,*), givptr(*), iwork(*), k(*), perm(ldgcol,*)
real(kind=BLASR8) :: c(*), d(*), difl(ldu,*), difr(ldu,*), e(*), givnum(ldu,*), poles(ldu,*), s(*), u(ldu,*), &
          vt(ldu,*), work(*), z(ldu,*)
call lb_dlasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z,poles,givptr,givcol,ldgcol,perm,givnum,c,s, &
               work,iwork,info)
end subroutine dlasda

subroutine dlasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)
use link_blas, only: lb_dlasdq
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, ldc, ldu, ldvt, n, ncc, ncvt, nru, sqre
real(kind=BLASR8) :: c(ldc,*), d(*), e(*), u(ldu,*), vt(ldvt,*), work(*)
call lb_dlasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)
end subroutine dlasdq

subroutine dlasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
use link_blas, only: lb_dlasdt
use Definitions, only: BLASInt
implicit none
integer(kind=BLASInt) :: lvl, msub, n, nd
integer(kind=BLASInt) :: inode(*), ndiml(*), ndimr(*)
call lb_dlasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
end subroutine dlasdt

subroutine dlaset(uplo,m,n,alpha,beta,a,lda)
use link_blas, only: lb_dlaset
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, m, n
real(kind=BLASR8) :: alpha, beta
real(kind=BLASR8) :: a(lda,*)
call lb_dlaset(uplo,m,n,alpha,beta,a,lda)
end subroutine dlaset

subroutine dlasq1(n,d,e,work,info)
use link_blas, only: lb_dlasq1
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: d(*), e(*), work(*)
call lb_dlasq1(n,d,e,work,info)
end subroutine dlasq1

subroutine dlasq2(n,z,info)
use link_blas, only: lb_dlasq2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: z(*)
call lb_dlasq2(n,z,info)
end subroutine dlasq2

subroutine dlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv,ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
use link_blas, only: lb_dlasq3
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: ieee
integer(kind=BLASInt) :: i0, iter, n0, ndiv, nfail, pp, ttype
real(kind=BLASR8) :: desig, dmin, dmin1, dmin2, dn, dn1, dn2, g, qmax, sigma, tau
real(kind=BLASR8) :: z(*)
call lb_dlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv,ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
end subroutine dlasq3

subroutine dlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype,g)
use link_blas, only: lb_dlasq4
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i0, n0, n0in, pp, ttype
real(kind=BLASR8) :: dmin, dmin1, dmin2, dn, dn1, dn2, g, tau
real(kind=BLASR8) :: z(*)
call lb_dlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype,g)
end subroutine dlasq4

subroutine dlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2,ieee,eps)
use link_blas, only: lb_dlasq5
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: ieee
integer(kind=BLASInt) :: i0, n0, pp
real(kind=BLASR8) :: dmin, dmin1, dmin2, dn, dnm1, dnm2, tau, sigma, eps
real(kind=BLASR8) :: z(*)
call lb_dlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2,ieee,eps)
end subroutine dlasq5

subroutine dlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
use link_blas, only: lb_dlasq6
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: i0, n0, pp
real(kind=BLASR8) :: dmin, dmin1, dmin2, dn, dnm1, dnm2
real(kind=BLASR8) :: z(*)
call lb_dlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
end subroutine dlasq6

subroutine dlasr(side,pivot,direct,m,n,c,s,a,lda)
use link_blas, only: lb_dlasr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, pivot, side
integer(kind=BLASInt) :: lda, m, n
real(kind=BLASR8) :: a(lda,*), c(*), s(*)
call lb_dlasr(side,pivot,direct,m,n,c,s,a,lda)
end subroutine dlasr

subroutine dlasrt(id,n,d,info)
use link_blas, only: lb_dlasrt
use Definitions, only: BLASInt, BLASR8
implicit none
character :: id
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: d(*)
call lb_dlasrt(id,n,d,info)
end subroutine dlasrt

subroutine dlassq(n,x,incx,scale,sumsq)
use link_blas, only: lb_dlassq
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: scale, sumsq
real(kind=BLASR8) :: x(*)
call lb_dlassq(n,x,incx,scale,sumsq)
end subroutine dlassq

subroutine dlasv2(f,g,h,ssmin,ssmax,snr,csr,snl,csl)
use link_blas, only: lb_dlasv2
use Definitions, only: BLASR8
implicit none
real(kind=BLASR8) :: csl, csr, f, g, h, snl, snr, ssmax, ssmin
call lb_dlasv2(f,g,h,ssmin,ssmax,snr,csr,snl,csl)
end subroutine dlasv2

subroutine dlaswp(n,a,lda,k1,k2,ipiv,incx)
use link_blas, only: lb_dlaswp
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, k1, k2, lda, n
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*)
call lb_dlaswp(n,a,lda,k1,k2,ipiv,incx)
end subroutine dlaswp

subroutine dlasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb,scale,x,ldx,xnorm,info)
use link_blas, only: lb_dlasy2
use Definitions, only: BLASInt, BLASR8
implicit none
logical :: ltranl, ltranr
integer(kind=BLASInt) :: info, isgn, ldb, ldtl, ldtr, ldx, n1, n2
real(kind=BLASR8) :: scale, xnorm
real(kind=BLASR8) :: b(ldb,*), tl(ldtl,*), tr(ldtr,*), x(ldx,*)
call lb_dlasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb,scale,x,ldx,xnorm,info)
end subroutine dlasy2

subroutine dlatrd(uplo,n,nb,a,lda,e,tau,w,ldw)
use link_blas, only: lb_dlatrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, ldw, n, nb
real(kind=BLASR8) :: a(lda,*), e(*), tau(*), w(ldw,*)
call lb_dlatrd(uplo,n,nb,a,lda,e,tau,w,ldw)
end subroutine dlatrd

subroutine dlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
use link_blas, only: lb_dlatrs
use Definitions, only: BLASInt, BLASR8
implicit none
character :: diag, normin, trans, uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: scale
real(kind=BLASR8) :: a(lda,*), cnorm(*), x(*)
call lb_dlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
end subroutine dlatrs

subroutine dopgtr(uplo,n,ap,tau,q,ldq,work,info)
use link_blas, only: lb_dopgtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, ldq, n
real(kind=BLASR8) :: ap(*), q(ldq,*), tau(*), work(*)
call lb_dopgtr(uplo,n,ap,tau,q,ldq,work,info)
end subroutine dopgtr

subroutine dopmtr(side,uplo,trans,m,n,ap,tau,c,ldc,work,info)
use link_blas, only: lb_dopmtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans, uplo
integer(kind=BLASInt) :: info, ldc, m, n
real(kind=BLASR8) :: ap(*), c(ldc,*), tau(*), work(*)
call lb_dopmtr(side,uplo,trans,m,n,ap,tau,c,ldc,work,info)
end subroutine dopmtr

subroutine dorg2l(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_dorg2l
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorg2l(m,n,k,a,lda,tau,work,info)
end subroutine dorg2l

subroutine dorg2r(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_dorg2r
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorg2r(m,n,k,a,lda,tau,work,info)
end subroutine dorg2r

subroutine dorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorgbr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: vect
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
end subroutine dorgbr

subroutine dorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorghr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: ihi, ilo, info, lda, lwork, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
end subroutine dorghr

subroutine dorgl2(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_dorgl2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorgl2(m,n,k,a,lda,tau,work,info)
end subroutine dorgl2

subroutine dorglq(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorglq
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorglq(m,n,k,a,lda,tau,work,lwork,info)
end subroutine dorglq

subroutine dorgql(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorgql
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorgql(m,n,k,a,lda,tau,work,lwork,info)
end subroutine dorgql

subroutine dorgqr(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorgqr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorgqr(m,n,k,a,lda,tau,work,lwork,info)
end subroutine dorgqr

subroutine dorgtr(uplo,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_dorgtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, lwork, n
real(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_dorgtr(uplo,n,a,lda,tau,work,lwork,info)
end subroutine dorgtr

subroutine dorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_dorm2l
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine dorm2l

subroutine dorm2r(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_dorm2r
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dorm2r(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine dorm2r

subroutine dormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormbr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans, vect
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormbr

subroutine dormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormhr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: ihi, ilo, info, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormhr

subroutine dorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_dorml2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine dorml2

subroutine dormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormlq
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormlq

subroutine dormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormql
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormql

subroutine dormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormqr

subroutine dormtr(side,uplo,trans,m,n,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_dormtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans, uplo
integer(kind=BLASInt) :: info, lda, ldc, lwork, m, n
real(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_dormtr(side,uplo,trans,m,n,a,lda,tau,c,ldc,work,lwork,info)
end subroutine dormtr

subroutine dposv(uplo,n,nrhs,a,lda,b,ldb,info)
use link_blas, only: lb_dposv
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dposv(uplo,n,nrhs,a,lda,b,ldb,info)
end subroutine dposv

subroutine dpotrf(uplo,n,a,lda,info)
use link_blas, only: lb_dpotrf
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*)
call lb_dpotrf(uplo,n,a,lda,info)
end subroutine dpotrf

recursive subroutine dpotrf2(uplo,n,a,lda,info)
use link_blas, only: lb_dpotrf2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*)
call lb_dpotrf2(uplo,n,a,lda,info)
end subroutine dpotrf2

subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
use link_blas, only: lb_dpotrs
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
end subroutine dpotrs

subroutine dpptrf(uplo,n,ap,info)
use link_blas, only: lb_dpptrf
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: ap(*)
call lb_dpptrf(uplo,n,ap,info)
end subroutine dpptrf

subroutine dpstf2(uplo,n,a,lda,piv,rank,tol,work,info)
use link_blas, only: lb_dpstf2
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: tol
integer(kind=BLASInt) :: info, lda, n, rank
character :: uplo
real(kind=BLASR8) :: a(lda,*), work(2*n)
integer(kind=BLASInt) :: piv(n)
call lb_dpstf2(uplo,n,a,lda,piv,rank,tol,work,info)
end subroutine dpstf2

subroutine dpstrf(uplo,n,a,lda,piv,rank,tol,work,info)
use link_blas, only: lb_dpstrf
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: tol
integer(kind=BLASInt) :: info, lda, n, rank
character :: uplo
real(kind=BLASR8) :: a(lda,*), work(2*n)
integer(kind=BLASInt) :: piv(n)
call lb_dpstrf(uplo,n,a,lda,piv,rank,tol,work,info)
end subroutine dpstrf

function droundup_lwork(lwork)
use link_blas, only: lb_droundup_lwork
use Definitions, only: BLASInt, BLASR8
implicit none
real(kind=BLASR8) :: droundup_lwork
integer(kind=BLASInt) :: lwork
droundup_lwork = lb_droundup_lwork(lwork)
end function droundup_lwork

subroutine drscl(n,sa,sx,incx)
use link_blas, only: lb_drscl
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: sa
real(kind=BLASR8) :: sx(*)
call lb_drscl(n,sa,sx,incx)
end subroutine drscl

subroutine dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
use link_blas, only: lb_dspev
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, ldz, n
real(kind=BLASR8) :: ap(*), w(*), work(*), z(ldz,*)
call lb_dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
end subroutine dspev

subroutine dspgst(itype,uplo,n,ap,bp,info)
use link_blas, only: lb_dspgst
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, itype, n
real(kind=BLASR8) :: ap(*), bp(*)
call lb_dspgst(itype,uplo,n,ap,bp,info)
end subroutine dspgst

subroutine dspgv(itype,jobz,uplo,n,ap,bp,w,z,ldz,work,info)
use link_blas, only: lb_dspgv
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, itype, ldz, n
real(kind=BLASR8) :: ap(*), bp(*), w(*), work(*), z(ldz,*)
call lb_dspgv(itype,jobz,uplo,n,ap,bp,w,z,ldz,work,info)
end subroutine dspgv

subroutine dsptrd(uplo,n,ap,d,e,tau,info)
use link_blas, only: lb_dsptrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: ap(*), d(*), e(*), tau(*)
call lb_dsptrd(uplo,n,ap,d,e,tau,info)
end subroutine dsptrd

subroutine dstebz(range,order,n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info)
use link_blas, only: lb_dstebz
use Definitions, only: BLASInt, BLASR8
implicit none
character :: order, range
integer(kind=BLASInt) :: il, info, iu, m, n, nsplit
real(kind=BLASR8) :: abstol, vl, vu
integer(kind=BLASInt) :: iblock(*), isplit(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), w(*), work(*)
call lb_dstebz(range,order,n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info)
end subroutine dstebz

subroutine dstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dstedc
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compz
integer(kind=BLASInt) :: info, ldz, liwork, lwork, n
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e(*), work(*), z(ldz,*)
call lb_dstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
end subroutine dstedc

subroutine dstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail,info)
use link_blas, only: lb_dstein
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, ldz, m, n
integer(kind=BLASInt) :: iblock(*), ifail(*), isplit(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), w(*), work(*), z(ldz,*)
call lb_dstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail,info)
end subroutine dstein

subroutine dstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc,isuppz,tryrac,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dstemr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, range
logical :: tryrac
integer(kind=BLASInt) :: il, info, iu, ldz, nzc, liwork, lwork, m, n
real(kind=BLASR8) :: vl, vu
integer(kind=BLASInt) :: isuppz(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), w(*), work(*)
real(kind=BLASR8) :: z(ldz,*)
call lb_dstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc,isuppz,tryrac,work,lwork,iwork,liwork,info)
end subroutine dstemr

subroutine dsteqr(compz,n,d,e,z,ldz,work,info)
use link_blas, only: lb_dsteqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compz
integer(kind=BLASInt) :: info, ldz, n
real(kind=BLASR8) :: d(*), e(*), work(*), z(ldz,*)
call lb_dsteqr(compz,n,d,e,z,ldz,work,info)
end subroutine dsteqr

subroutine dsterf(n,d,e,info)
use link_blas, only: lb_dsterf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: d(*), e(*)
call lb_dsterf(n,d,e,info)
end subroutine dsterf

subroutine dstevr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dstevr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, range
integer(kind=BLASInt) :: il, info, iu, ldz, liwork, lwork, m, n
real(kind=BLASR8) :: abstol, vl, vu
integer(kind=BLASInt) :: isuppz(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), w(*), work(*), z(ldz,*)
call lb_dstevr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
end subroutine dstevr

subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
use link_blas, only: lb_dsyev
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, lda, lwork, n
real(kind=BLASR8) :: a(lda,*), w(*), work(*)
call lb_dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
end subroutine dsyev

subroutine dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dsyevd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, lda, liwork, lwork, n
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: a(lda,*), w(*), work(*)
call lb_dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
end subroutine dsyevd

subroutine dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dsyevr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, range, uplo
integer(kind=BLASInt) :: il, info, iu, lda, ldz, liwork, lwork, m, n
real(kind=BLASR8) :: abstol, vl, vu
integer(kind=BLASInt) :: isuppz(*), iwork(*)
real(kind=BLASR8) :: a(lda,*), w(*), work(*), z(ldz,*)
call lb_dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
end subroutine dsyevr

subroutine dsygs2(itype,uplo,n,a,lda,b,ldb,info)
use link_blas, only: lb_dsygs2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, itype, lda, ldb, n
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dsygs2(itype,uplo,n,a,lda,b,ldb,info)
end subroutine dsygs2

subroutine dsygst(itype,uplo,n,a,lda,b,ldb,info)
use link_blas, only: lb_dsygst
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, itype, lda, ldb, n
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dsygst(itype,uplo,n,a,lda,b,ldb,info)
end subroutine dsygst

subroutine dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
use link_blas, only: lb_dsygv
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, itype, lda, ldb, lwork, n
real(kind=BLASR8) :: a(lda,*), b(ldb,*), w(*), work(*)
call lb_dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
end subroutine dsygv

subroutine dsytd2(uplo,n,a,lda,d,e,tau,info)
use link_blas, only: lb_dsytd2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*), d(*), e(*), tau(*)
call lb_dsytd2(uplo,n,a,lda,d,e,tau,info)
end subroutine dsytd2

subroutine dsytrd(uplo,n,a,lda,d,e,tau,work,lwork,info)
use link_blas, only: lb_dsytrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, lwork, n
real(kind=BLASR8) :: a(lda,*), d(*), e(*), tau(*), work(*)
call lb_dsytrd(uplo,n,a,lda,d,e,tau,work,lwork,info)
end subroutine dsytrd

subroutine dtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m,work,lwork,info)
use link_blas, only: lb_dtrevc3
use Definitions, only: BLASInt, BLASR8
implicit none
character :: howmny, side
integer(kind=BLASInt) :: info, ldt, ldvl, ldvr, lwork, m, mm, n
logical :: select(*)
real(kind=BLASR8) :: t(ldt,*), vl(ldvl,*), vr(ldvr,*), work(*)
call lb_dtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m,work,lwork,info)
end subroutine dtrevc3

subroutine dtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
use link_blas, only: lb_dtrexc
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compq
integer(kind=BLASInt) :: ifst, ilst, info, ldq, ldt, n
real(kind=BLASR8) :: q(ldq,*), t(ldt,*), work(*)
call lb_dtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
end subroutine dtrexc

subroutine dtrsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work,lwork,iwork,liwork,info)
use link_blas, only: lb_dtrsen
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compq, job
integer(kind=BLASInt) :: info, ldq, ldt, liwork, lwork, m, n
real(kind=BLASR8) :: s, sep
logical :: select(*)
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: q(ldq,*), t(ldt,*), wi(*), work(*), wr(*)
call lb_dtrsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work,lwork,iwork,liwork,info)
end subroutine dtrsen

subroutine dtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
use link_blas, only: lb_dtrsyl
use Definitions, only: BLASInt, BLASR8
implicit none
character :: trana, tranb
integer(kind=BLASInt) :: info, isgn, lda, ldb, ldc, m, n
real(kind=BLASR8) :: scale
real(kind=BLASR8) :: a(lda,*), b(ldb,*), c(ldc,*)
call lb_dtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
end subroutine dtrsyl

subroutine dtrti2(uplo,diag,n,a,lda,info)
use link_blas, only: lb_dtrti2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: diag, uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*)
call lb_dtrti2(uplo,diag,n,a,lda,info)
end subroutine dtrti2

subroutine dtrtri(uplo,diag,n,a,lda,info)
use link_blas, only: lb_dtrtri
use Definitions, only: BLASInt, BLASR8
implicit none
character :: diag, uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*)
call lb_dtrtri(uplo,diag,n,a,lda,info)
end subroutine dtrtri

subroutine dtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
use link_blas, only: lb_dtrtrs
use Definitions, only: BLASInt, BLASR8
implicit none
character :: diag, trans, uplo
integer(kind=BLASInt) :: info, lda, ldb, n, nrhs
real(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_dtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
end subroutine dtrtrs

function ieeeck(ispec,zero,one)
use link_blas, only: lb_ieeeck
use Definitions, only: BLASInt, BLASR4
implicit none
integer(kind=BLASInt) :: ispec
real(kind=BLASR4) :: one, zero
integer(kind=BLASInt) :: ieeeck
ieeeck = lb_ieeeck(ispec,zero,one)
end function ieeeck

function iladlc(m,n,a,lda)
use link_blas, only: lb_iladlc
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: m, n, lda
real(kind=BLASR8) :: a(lda,*)
integer(kind=BLASInt) :: iladlc
iladlc = lb_iladlc(m,n,a,lda)
end function iladlc

function iladlr(m,n,a,lda)
use link_blas, only: lb_iladlr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: m, n, lda
real(kind=BLASR8) :: a(lda,*)
integer(kind=BLASInt) :: iladlr
iladlr = lb_iladlr(m,n,a,lda)
end function iladlr

function ilaenv(ispec,name,opts,n1,n2,n3,n4)
use link_blas, only: lb_ilaenv
use Definitions, only: BLASInt
implicit none
character(len=*) :: name, opts
integer(kind=BLASInt) :: ispec, n1, n2, n3, n4
integer(kind=BLASInt) :: ilaenv
ilaenv = lb_ilaenv(ispec,name,opts,n1,n2,n3,n4)
end function ilaenv

function ilazlc(m,n,a,lda)
use link_blas, only: lb_ilazlc
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: m, n, lda
complex(kind=BLASR8) :: a(lda,*)
integer(kind=BLASInt) :: ilazlc
ilazlc = lb_ilazlc(m,n,a,lda)
end function ilazlc

function ilazlr(m,n,a,lda)
use link_blas, only: lb_ilazlr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: m, n, lda
complex(kind=BLASR8) :: a(lda,*)
integer(kind=BLASInt) :: ilazlr
ilazlr = lb_ilazlr(m,n,a,lda)
end function ilazlr

function iparam2stage(ispec,name,opts,ni,nbi,ibi,nxi)
use link_blas, only: lb_iparam2stage
use Definitions, only: BLASInt
implicit none
character(len=*) :: name, opts
integer(kind=BLASInt) :: ispec, ni, nbi, ibi, nxi
integer(kind=BLASInt) :: iparam2stage
iparam2stage = lb_iparam2stage(ispec,name,opts,ni,nbi,ibi,nxi)
end function iparam2stage

function iparmq(ispec,name,opts,n,ilo,ihi,lwork)
use link_blas, only: lb_iparmq
use Definitions, only: BLASInt
implicit none
integer(kind=BLASInt) :: ihi, ilo, ispec, lwork, n
character(len=*) :: name, opts
integer(kind=BLASInt) :: iparmq
iparmq = lb_iparmq(ispec,name,opts,n,ilo,ihi,lwork)
end function iparmq

subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
use link_blas, only: lb_zheev
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, lda, lwork, n
real(kind=BLASR8) :: rwork(*), w(*)
complex(kind=BLASR8) :: a(lda,*), work(*)
call lb_zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
end subroutine zheev

subroutine zhetd2(uplo,n,a,lda,d,e,tau,info)
use link_blas, only: lb_zhetd2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: a(lda,*), tau(*)
call lb_zhetd2(uplo,n,a,lda,d,e,tau,info)
end subroutine zhetd2

subroutine zhetrd(uplo,n,a,lda,d,e,tau,work,lwork,info)
use link_blas, only: lb_zhetrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, lwork, n
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zhetrd(uplo,n,a,lda,d,e,tau,work,lwork,info)
end subroutine zhetrd

subroutine zhpev(jobz,uplo,n,ap,w,z,ldz,work,rwork,info)
use link_blas, only: lb_zhpev
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobz, uplo
integer(kind=BLASInt) :: info, ldz, n
real(kind=BLASR8) :: rwork(*), w(*)
complex(kind=BLASR8) :: ap(*), work(*), z(ldz,*)
call lb_zhpev(jobz,uplo,n,ap,w,z,ldz,work,rwork,info)
end subroutine zhpev

subroutine zhptrd(uplo,n,ap,d,e,tau,info)
use link_blas, only: lb_zhptrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, n
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: ap(*), tau(*)
call lb_zhptrd(uplo,n,ap,d,e,tau,info)
end subroutine zhptrd

subroutine zlacgv(n,x,incx)
use link_blas, only: lb_zlacgv
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
complex(kind=BLASR8) :: x(*)
call lb_zlacgv(n,x,incx)
end subroutine zlacgv

function zladiv(x,y)
use link_blas, only: lb_zladiv
use Definitions, only: BLASR8
implicit none
complex(kind=BLASR8) :: x, y
complex(kind=BLASR8) :: zladiv
zladiv = lb_zladiv(x,y)
end function zladiv

function zlange(norm,m,n,a,lda,work)
use link_blas, only: lb_zlange
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm
integer(kind=BLASInt) :: lda, m, n
real(kind=BLASR8) :: work(*)
complex(kind=BLASR8) :: a(lda,*)
real(kind=BLASR8) :: zlange
zlange = lb_zlange(norm,m,n,a,lda,work)
end function zlange

function zlanhe(norm,uplo,n,a,lda,work)
use link_blas, only: lb_zlanhe
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm, uplo
integer(kind=BLASInt) :: lda, n
real(kind=BLASR8) :: work(*)
complex(kind=BLASR8) :: a(lda,*)
real(kind=BLASR8) :: zlanhe
zlanhe = lb_zlanhe(norm,uplo,n,a,lda,work)
end function zlanhe

function zlanhp(norm,uplo,n,ap,work)
use link_blas, only: lb_zlanhp
use Definitions, only: BLASInt, BLASR8
implicit none
character :: norm, uplo
integer(kind=BLASInt) :: n
real(kind=BLASR8) :: work(*)
complex(kind=BLASR8) :: ap(*)
real(kind=BLASR8) :: zlanhp
zlanhp = lb_zlanhp(norm,uplo,n,ap,work)
end function zlanhp

subroutine zlarf(side,m,n,v,incv,tau,c,ldc,work)
use link_blas, only: lb_zlarf
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side
integer(kind=BLASInt) :: incv, ldc, m, n
complex(kind=BLASR8) :: tau
complex(kind=BLASR8) :: c(ldc,*), v(*), work(*)
call lb_zlarf(side,m,n,v,incv,tau,c,ldc,work)
end subroutine zlarf

subroutine zlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc,work,ldwork)
use link_blas, only: lb_zlarfb
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, side, storev, trans
integer(kind=BLASInt) :: k, ldc, ldt, ldv, ldwork, m, n
complex(kind=BLASR8) :: c(ldc,*), t(ldt,*), v(ldv,*), work(ldwork,*)
call lb_zlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc,work,ldwork)
end subroutine zlarfb

subroutine zlarfg(n,alpha,x,incx,tau)
use link_blas, only: lb_zlarfg
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
complex(kind=BLASR8) :: alpha, tau
complex(kind=BLASR8) :: x(*)
call lb_zlarfg(n,alpha,x,incx,tau)
end subroutine zlarfg

subroutine zlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
use link_blas, only: lb_zlarft
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, storev
integer(kind=BLASInt) :: k, ldt, ldv, n
complex(kind=BLASR8) :: t(ldt,*), tau(*), v(ldv,*)
call lb_zlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
end subroutine zlarft

subroutine zlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
use link_blas, only: lb_zlascl
use Definitions, only: BLASInt, BLASR8
implicit none
character :: type
integer(kind=BLASInt) :: info, kl, ku, lda, m, n
real(kind=BLASR8) :: cfrom, cto
complex(kind=BLASR8) :: a(lda,*)
call lb_zlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
end subroutine zlascl

subroutine zlaset(uplo,m,n,alpha,beta,a,lda)
use link_blas, only: lb_zlaset
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, m, n
complex(kind=BLASR8) :: alpha, beta
complex(kind=BLASR8) :: a(lda,*)
call lb_zlaset(uplo,m,n,alpha,beta,a,lda)
end subroutine zlaset

subroutine zlasr(side,pivot,direct,m,n,c,s,a,lda)
use link_blas, only: lb_zlasr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: direct, pivot, side
integer(kind=BLASInt) :: lda, m, n
real(kind=BLASR8) :: c(*), s(*)
complex(kind=BLASR8) :: a(lda,*)
call lb_zlasr(side,pivot,direct,m,n,c,s,a,lda)
end subroutine zlasr

subroutine zlassq(n,x,incx,scale,sumsq)
use link_blas, only: lb_zlassq
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: incx, n
real(kind=BLASR8) :: scale, sumsq
complex(kind=BLASR8) :: x(*)
call lb_zlassq(n,x,incx,scale,sumsq)
end subroutine zlassq

subroutine zlatrd(uplo,n,nb,a,lda,e,tau,w,ldw)
use link_blas, only: lb_zlatrd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, ldw, n, nb
real(kind=BLASR8) :: e(*)
complex(kind=BLASR8) :: a(lda,*), tau(*), w(ldw,*)
call lb_zlatrd(uplo,n,nb,a,lda,e,tau,w,ldw)
end subroutine zlatrd

subroutine zsteqr(compz,n,d,e,z,ldz,work,info)
use link_blas, only: lb_zsteqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compz
integer(kind=BLASInt) :: info, ldz, n
real(kind=BLASR8) :: d(*), e(*), work(*)
complex(kind=BLASR8) :: z(ldz,*)
call lb_zsteqr(compz,n,d,e,z,ldz,work,info)
end subroutine zsteqr

subroutine zung2l(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_zung2l
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zung2l(m,n,k,a,lda,tau,work,info)
end subroutine zung2l

subroutine zung2r(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_zung2r
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zung2r(m,n,k,a,lda,tau,work,info)
end subroutine zung2r

subroutine zungql(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zungql
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zungql(m,n,k,a,lda,tau,work,lwork,info)
end subroutine zungql

subroutine zungqr(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zungqr
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zungqr(m,n,k,a,lda,tau,work,lwork,info)
end subroutine zungqr

subroutine zungtr(uplo,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zungtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, lwork, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zungtr(uplo,n,a,lda,tau,work,lwork,info)
end subroutine zungtr

subroutine zupgtr(uplo,n,ap,tau,q,ldq,work,info)
use link_blas, only: lb_zupgtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, ldq, n
complex(kind=BLASR8) :: ap(*), q(ldq,*), tau(*), work(*)
call lb_zupgtr(uplo,n,ap,tau,q,ldq,work,info)
end subroutine zupgtr

subroutine zlacpy(uplo,m,n,a,lda,b,ldb)
use link_blas, only: lb_zlacpy
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: lda, ldb, m, n
complex(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_zlacpy(uplo,m,n,a,lda,b,ldb)
end subroutine zlacpy

subroutine zbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork,info)
use link_blas, only: lb_zbdsqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, ldc, ldu, ldvt, n, ncc, ncvt, nru
real(kind=BLASR8) :: d(*), e(*), rwork(*)
complex(kind=BLASR8) :: c(ldc,*), u(ldu,*), vt(ldvt,*)
call lb_zbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork,info)
end subroutine zbdsqr

subroutine zgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
use link_blas, only: lb_zgebd2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: a(lda,*), taup(*), tauq(*), work(*)
call lb_zgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
end subroutine zgebd2

subroutine zgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
use link_blas, only: lb_zgebrd
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: a(lda,*), taup(*), tauq(*), work(*)
call lb_zgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
end subroutine zgebrd

subroutine zgelq2(m,n,a,lda,tau,work,info)
use link_blas, only: lb_zgelq2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zgelq2(m,n,a,lda,tau,work,info)
end subroutine zgelq2

subroutine zgelqf(m,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zgelqf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zgelqf(m,n,a,lda,tau,work,lwork,info)
end subroutine zgelqf

subroutine zgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
use link_blas, only: lb_zgels
use Definitions, only: BLASInt, BLASR8
implicit none
character :: trans
integer(kind=BLASInt) :: m, n, nrhs, lda, ldb, lwork, info
complex(kind=BLASR8) :: a(lda,*), b(ldb,*), work(*)
call lb_zgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
end subroutine zgels

subroutine zgeqr2(m,n,a,lda,tau,work,info)
use link_blas, only: lb_zgeqr2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zgeqr2(m,n,a,lda,tau,work,info)
end subroutine zgeqr2

subroutine zgeqrf(m,n,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zgeqrf
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zgeqrf(m,n,a,lda,tau,work,lwork,info)
end subroutine zgeqrf

subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
use link_blas, only: lb_zgesvd
use Definitions, only: BLASInt, BLASR8
implicit none
character :: jobu, jobvt
integer(kind=BLASInt) :: info, lda, ldu, ldvt, lwork, m, n
real(kind=BLASR8) :: rwork(*), s(*)
complex(kind=BLASR8) :: a(lda,*), u(ldu,*), vt(ldvt,*), work(*)
call lb_zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
end subroutine zgesvd

subroutine zlabrd(m,n,nb,a,lda,d,e,tauq,taup,x,ldx,y,ldy)
use link_blas, only: lb_zlabrd
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: lda, ldx, ldy, m, n, nb
real(kind=BLASR8) :: d(*), e(*)
complex(kind=BLASR8) :: a(lda,*), taup(*), tauq(*), x(ldx,*), y(ldy,*)
call lb_zlabrd(m,n,nb,a,lda,d,e,tauq,taup,x,ldx,y,ldy)
end subroutine zlabrd

subroutine zlacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
use link_blas, only: lb_zlacrm
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: lda, ldb, ldc, m, n
real(kind=BLASR8) :: b(ldb,*), rwork(*)
complex(kind=BLASR8) :: a(lda,*), c(ldc,*)
call lb_zlacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
end subroutine zlacrm

subroutine zlaed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
use link_blas, only: lb_zlaed0
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, ldq, ldqs, n, qsiz
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e(*), rwork(*)
complex(kind=BLASR8) :: q(ldq,*), qstore(ldqs,*)
call lb_zlaed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
end subroutine zlaed0

subroutine zlaed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq,qstore,qptr,prmptr,perm,givptr,givcol, &
                  givnum,work,rwork,iwork,info)
use link_blas, only: lb_zlaed7
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: curlvl, curpbm, cutpnt, info, ldq, n, qsiz, tlvls
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: givcol(2,*), givptr(*), indxq(*), iwork(*), perm(*), prmptr(*), qptr(*)
real(kind=BLASR8) :: d(*), givnum(2,*), qstore(*), rwork(*)
complex(kind=BLASR8) :: q(ldq,*), work(*)
call lb_zlaed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq,qstore,qptr,prmptr,perm,givptr,givcol, &
               givnum,work,rwork,iwork,info)
end subroutine zlaed7

subroutine zlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w,indxp,indx,indxq,perm,givptr,givcol,givnum, &
                  info)
use link_blas, only: lb_zlaed8
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: cutpnt, givptr, info, k, ldq, ldq2, n, qsiz
real(kind=BLASR8) :: rho
integer(kind=BLASInt) :: givcol(2,*), indx(*), indxp(*), indxq(*), perm(*)
real(kind=BLASR8) :: d(*), dlamda(*), givnum(2,*), w(*), z(*)
complex(kind=BLASR8) :: q(ldq,*), q2(ldq2,*)
call lb_zlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w,indxp,indx,indxq,perm,givptr,givcol,givnum, &
               info)
end subroutine zlaed8

subroutine zstedc(compz,n,d,e,z,ldz,work,lwork,rwork,lrwork,iwork,liwork,info)
use link_blas, only: lb_zstedc
use Definitions, only: BLASInt, BLASR8
implicit none
character :: compz
integer(kind=BLASInt) :: info, ldz, liwork, lrwork, lwork, n
integer(kind=BLASInt) :: iwork(*)
real(kind=BLASR8) :: d(*), e(*), rwork(*)
complex(kind=BLASR8) :: work(*), z(ldz,*)
call lb_zstedc(compz,n,d,e,z,ldz,work,lwork,rwork,lrwork,iwork,liwork,info)
end subroutine zstedc

subroutine zstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail,info)
use link_blas, only: lb_zstein
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, ldz, m, n
integer(kind=BLASInt) :: iblock(*), ifail(*), isplit(*), iwork(*)
real(kind=BLASR8) :: d(*), e(*), w(*), work(*)
complex(kind=BLASR8) :: z(ldz,*)
call lb_zstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail,info)
end subroutine zstein

subroutine ztrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
use link_blas, only: lb_ztrtrs
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo, trans, diag
integer(kind=BLASInt) :: n, nrhs, lda, ldb, info
complex(kind=BLASR8) :: a(lda,*), b(ldb,*)
call lb_ztrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
end subroutine ztrtrs

subroutine zungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zungbr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: vect
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
end subroutine zungbr

subroutine zungl2(m,n,k,a,lda,tau,work,info)
use link_blas, only: lb_zungl2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zungl2(m,n,k,a,lda,tau,work,info)
end subroutine zungl2

subroutine zunglq(m,n,k,a,lda,tau,work,lwork,info)
use link_blas, only: lb_zunglq
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, k, lda, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), tau(*), work(*)
call lb_zunglq(m,n,k,a,lda,tau,work,lwork,info)
end subroutine zunglq

subroutine zunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_zunm2l
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine zunm2l

subroutine zunm2r(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_zunm2r
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunm2r(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine zunm2r

subroutine zunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_zunmbr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans, vect
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine zunmbr

subroutine zunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
use link_blas, only: lb_zunml2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine zunml2

subroutine zunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_zunmlq
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine zunmlq

subroutine zunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_zunmql
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine zunmql

subroutine zunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_zunmqr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans
integer(kind=BLASInt) :: info, k, lda, ldc, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)
end subroutine zunmqr

subroutine zunmtr(side,uplo,trans,m,n,a,lda,tau,c,ldc,work,lwork,info)
use link_blas, only: lb_zunmtr
use Definitions, only: BLASInt, BLASR8
implicit none
character :: side, trans, uplo
integer(kind=BLASInt) :: info, lda, ldc, lwork, m, n
complex(kind=BLASR8) :: a(lda,*), c(ldc,*), tau(*), work(*)
call lb_zunmtr(side,uplo,trans,m,n,a,lda,tau,c,ldc,work,lwork,info)
end subroutine zunmtr
