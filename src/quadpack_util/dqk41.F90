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

subroutine dqk41(f,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk41
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  41-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real*8 version
!
!           parameters
!            on entry
!              f      - real*8
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real*8
!                       lower limit of integration
!
!              b      - real*8
!                       upper limit of integration
!
!            on return
!              reslt  - real*8
!                       approximation to the integral i
!                       result is computed by applying the 41-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 20-point gauss rule (resg).
!
!              abserr - real*8
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-reslt)
!
!              resabs - real*8
!                       approximation to the integral j
!
!              resasc - real*8
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk41

use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
interface
  function f(x)
    import :: wp
    real(kind=wp) :: f
    real(kind=wp), intent(in) :: x
  end function f
end interface
real(kind=wp), intent(in) :: a, b
real(kind=wp), intent(out) :: reslt, abserr, resabs, resasc
real(kind=wp) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(20), fv2(20), hlgth, resg, resk, reskh, uflow
integer(kind=iwp) :: j, jtw, jtwm1
real(kind=wp), external :: d1mach

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
! xgk    - abscissae of the 41-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 20-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 20-point gauss rule
!
! wgk    - weights of the 41-point kronrod rule
!
! wg     - weights of the 20-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

real(kind=wp), parameter :: wg(10) = [0.017614007139152118311861962351853_wp, &
                                      0.040601429800386941331039952274932_wp, &
                                      0.062672048334109063569506535187042_wp, &
                                      0.083276741576704748724758143222046_wp, &
                                      0.101930119817240435036750135480350_wp, &
                                      0.118194531961518417312377377711382_wp, &
                                      0.131688638449176626898494499748163_wp, &
                                      0.142096109318382051329298325067165_wp, &
                                      0.149172986472603746787828737001969_wp, &
                                      0.152753387130725850698084331955098_wp]
real(kind=wp), parameter :: xgk(21) = [0.998859031588277663838315576545863_wp, &
                                       0.993128599185094924786122388471320_wp, &
                                       0.981507877450250259193342994720217_wp, &
                                       0.963971927277913791267666131197277_wp, &
                                       0.940822633831754753519982722212443_wp, &
                                       0.912234428251325905867752441203298_wp, &
                                       0.878276811252281976077442995113078_wp, &
                                       0.839116971822218823394529061701521_wp, &
                                       0.795041428837551198350638833272788_wp, &
                                       0.746331906460150792614305070355642_wp, &
                                       0.693237656334751384805490711845932_wp, &
                                       0.636053680726515025452836696226286_wp, &
                                       0.575140446819710315342946036586425_wp, &
                                       0.510867001950827098004364050955251_wp, &
                                       0.443593175238725103199992213492640_wp, &
                                       0.373706088715419560672548177024927_wp, &
                                       0.301627868114913004320555356858592_wp, &
                                       0.227785851141645078080496195368575_wp, &
                                       0.152605465240922675505220241022678_wp, &
                                       0.076526521133497333754640409398838_wp, &
                                       0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(21) = [0.003073583718520531501218293246031_wp, &
                                       0.008600269855642942198661787950102_wp, &
                                       0.014626169256971252983787960308868_wp, &
                                       0.020388373461266523598010231432755_wp, &
                                       0.025882133604951158834505067096153_wp, &
                                       0.031287306777032798958543119323801_wp, &
                                       0.036600169758200798030557240707211_wp, &
                                       0.041668873327973686263788305936895_wp, &
                                       0.046434821867497674720231880926108_wp, &
                                       0.050944573923728691932707670050345_wp, &
                                       0.055195105348285994744832372419777_wp, &
                                       0.059111400880639572374967220648594_wp, &
                                       0.062653237554781168025870122174255_wp, &
                                       0.065834597133618422111563556969398_wp, &
                                       0.068648672928521619345623411885368_wp, &
                                       0.071054423553444068305790361723210_wp, &
                                       0.073030690332786667495189417658913_wp, &
                                       0.074582875400499188986581418362488_wp, &
                                       0.075704497684556674659542775376617_wp, &
                                       0.076377867672080736705502835038061_wp, &
                                       0.076600711917999656445049901530102_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc   - abscissa
! fval*  - function value
! resg   - result of the 20-point gauss formula
! resk   - result of the 41-point kronrod formula
! reskh  - approximation to the mean value of f over (a,b),
!          i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk41

epmach = d1mach(4)
uflow = d1mach(1)

centr = Half*(a+b)
hlgth = Half*(b-a)
dhlgth = abs(hlgth)

! compute the 41-point kronrod approximation to
! the integral, and estimate the absolute error.

fc = f(centr)
resg = Zero
resk = fc*wgk(21)
resabs = abs(resk)
do j=1,10
  jtw = j*2
  absc = hlgth*xgk(jtw)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1+fval2
  resg = resg+wg(j)*fsum
  resk = resk+wgk(jtw)*fsum
  resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
end do
do j=1,10
  jtwm1 = j*2-1
  absc = hlgth*xgk(jtwm1)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1+fval2
  resk = resk+wgk(jtwm1)*fsum
  resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
end do
reskh = resk*Half
resasc = wgk(21)*abs(fc-reskh)
do j=1,20
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*dhlgth
resasc = resasc*dhlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk41
