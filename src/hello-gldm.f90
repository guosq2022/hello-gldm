module hello_gldm
    implicit none
    private

    public :: say_hello, gldm

    real(8) :: temp = 0., paray = 0.02
    integer :: nr1 = 1, nmax = 800, npas = 1, n1cor = 101
    integer, parameter :: opt1 = 1 ! to write s, sp, c, es ec, en, et
    integer, parameter :: opt2 = 0 ! to give energy relative to the infinity (fusion barrier).
    integer, parameter :: opt3 = 0 ! to give energy relative to the ground state (fission barrier).
    integer, parameter :: opt4 = 1 ! to determine the radii from the total system, otherwise, from each nucleus.
    integer, parameter :: opt5 = 0 ! to calculate fusion cross sections using GLDM.
    integer, parameter :: opt6 = 1 ! to give ra, rb, tl, hw, elmax (characteristics of fusion barrier).
    integer, parameter :: opt7 = 1 ! to take into account ellipsoidal deformation for fusion.
    integer, parameter :: opt8 = 1 ! to take into account the shell effects
    integer, parameter :: opt10 = 0 ! to calculate fusion cross section using WONG FORMULA.
    integer, parameter :: opt11 = 0
    integer, parameter :: opt12 = 0

contains

    subroutine say_hello
        print *, "Hello, hello-gldm!"
        print *, " test github working!"
    end subroutine say_hello

    ! call gldm
    subroutine gldm(frag_a1, frag_z1, frag_a2, frag_z2, qexp_in, t12_out)
        real(8), intent(in) :: frag_a1, frag_z1, frag_a2, frag_z2, qexp_in
        real(8), intent(out) :: t12_out

        real(8) :: a0, a1, a2, as, avol, alpha, a1unt, a1det, a2unt, a2det, &
        & aunti, adeti, amo, aco, ac2, ac4, ac6, arc, amort, akr, ak2, aak, &
        & a, axxo, actn, auxse
        real(8) :: beta, beta1, beta2, beta3, bdif, bs, bc, b
        real(8) :: coefq, cosep, coec, c, c2, c3, c4, c5, c7, c8, c9, &
        & c2s2, cp2, cp, cp3, cp4, cpac, c2ac, csd, cpt12, cpzv, czv, &
        & cp4z, c4z, cos2x, couch, coooo, co2oo, cosoi, cos2, css, ccc, &
        & cpc, coooi, cor, cor1, cor2, coefk, coefd, cefd2, croi, coeff1, &
        & concd, conca, chera, cherb
        real(8) :: depi, difvol, difes, difec, disco, denbc, d5, d2, d, &
        & dx, devr1, devr2, de2, deooi, d9, d1, dvj, di4m1, dpent, deshw
        real(8) :: eso, ess, eco, ecinf, evoinf, evolsp, eo, eninf, eps, es, &
        & ecart, ec, en, e, etota, erot, ecer, ebarr, ebafi, econt, eref, eee, &
        & elmax, evme, emaxl
        real(8) :: fs, f, fl, fr, fm, fm1, fo, fp1, freq
        real(8) :: h, hh, hv, hy, h1, h2, h1ph2, hw, hww
        real(8) :: llll
        real(8) :: ooo, oi, ooooi
        real(8) :: pas, ph1, ph2, ph3, phs1, phs2, phs3, pisur2, pasin, &
        & penta, pentb, paeb, pp12, prodd, pila2, prod
        real(8) :: qexp, qvalu, qupi1, qupi2, qupig, qreac, quadr, q2
        real(8) :: r1, r2, ro, rkas, raya1, raya2, rfiss, rayon, rayo2, rayo3, &
        & r1pr2, rac, rat12, rcent, r13, r23, r2x, rrx, rrr, rcc, r2j, rjj, &
        & rayo5, rbarr, redmas, rol, rkl, rra, rrb, ra
        real(8) :: sep, s, s2, s4, sprim, sp2, sp4, sqzv, sqcp, sqc, &
        & sq1, sq2, sqp1, sqp2, ssur1, ssup1, sinox, sin2x, shell,   &
        & s2m1, sp2m1, siooo, sinoi, sin2, siooi, s6, step, secte, sect
        real(8) :: t12, tcrit, tol, tolf, tetax, tzhw, tvar1, test, tvar2, ttl
        real(8) :: r2k, rkk, r1t, r2t, rout, routth, rb, rrol
        real(8) :: uti, u
        real(8) :: vo1, vpot1, vpot2, vartl
        real(8) :: wx, w1, w2, w3, ws1, ws2, ws3, wh
        real(8) :: x, xuu, xli, xri, xl1, xr, xm, xcent, xl, xu, &
        & xxmoq, xmom2, xmo1, xmo2, xlmom, xlamb, xlog10
        real(8) :: y
        real(8) :: z0, z1, z2, z1i, z1i2, z2i, z2i2, zi, zii, z1z2, z2sua, z1sa1, &
        & z2sa2, z4, za, zv, zv2, zv4, zz1, zz2, zhw
        real(8), parameter :: pi = 3.1415926535
        real(8), dimension(2000) :: rai, etot, rt, erego, erefi, deri, et1, zz, tl
        real(8), dimension(100) :: rintx, devr, y1, vo, voo, bcc
        integer, dimension(10) :: nmin12
        integer :: nsec1, nsec2, nh1, nh2, nr1p1, nm, nmax1, i, j, k, iend, ier, m, &
        & npasin, intrx, jj, nema, nem, jr, ir, kr, jrt, krm1, jd, &
        & ncont, nminim, n111, n222, ij, nabsmi, ji, ikk, nturn, ii, l, ll, jm, ijk, &
        & ie, ll5, lr, jrm1, lmax, jol, lex, irai, l2, ndim, lrr, lvl, l7, jjol, &
        & nndim, le, l1

        open (10, file='result.dat', status='replace')

        a2 = frag_a2
        z2 = frag_z2
        z1 = frag_z1
        a1 = frag_a1
        qexp = qexp_in
        llll = 10.0
        xlmom = 10.0
        write (*, *) 'a1, z1, a2, z2, qexp:', a1, z1, a2, z2, qexp


        ! old main gldm
        a0 = a1 + a2
        z0 = z1 + z2

        tcrit = 17.
        qvalu = 0.0

        ! NSEC -7 to -7+NSEC1 STEP NSEC2 to have energy for fusiopn cross sections
        nsec1 = 160
        nsec2 = 8
        nh1 = 25
        nh2 = 29
        hv = 3.141470/(1.*(nh1 - 1)*1.)
        hy = 3.141290/(1.*(nh2 - 1)*1.)
        as = 17.9439*(1.+1.5*temp/tcrit)*((1.-temp/tcrit)**1.5)
        avol = -15.4941*(1.+0.00337*temp*temp)
        ro = 1.2249
        rkas = 2.6
        uti = 1./3.
        nr1p1 = nr1 + 1
        nm = nr1 + 1
        nmax1 = nmax - 1
        depi = 6.2831853
        h = 0.1
        hh = h*h
        alpha = (a1 - a2)/(a1 + a2)
        a1unt = a1**.3333333333
        a1det = a1**.6666666666
        a2unt = a2**.3333333333
        a2det = a2**.6666666666
        aunti = a0**.3333333333
        adeti = a0**.6666666666
        z1i = (a1 - 2.*z1)/a1
        z1i2 = z1i*z1i
        z2i = (a2 - 2.*z2)/a2
        z2i2 = z2i*z2i
        zi = (a0 - 2.*z0)/a0
        zii = zi*zi
        raya1 = 1.28*a1unt - 0.76 + (0.8/a1unt)
        raya2 = 1.28*a2unt - 0.76 + (0.8/a2unt)
        if (opt4 == 0) then
            beta1 = raya2/raya1
            beta3 = beta1*beta1*beta1
            rayon = 1.28*aunti - 0.76 + (0.8/aunti)
            rayon = rayon*(1.0 + 0.0007*temp*temp)
            raya1 = rayon/(1.0 + beta3)**0.33333333
            raya2 = beta1*raya1
        else
            rayon = (raya1**3.0 + raya2**3.0)**uti
            rayon = rayon*(1.0 + 0.0007*temp*temp)
            raya1 = raya1*(1.0 + 0.0007*temp*temp)
            raya2 = raya2*(1.0 + 0.0007*temp*temp)
        end if
        rayo2 = rayon*rayon

        ! Coefficient to express the quadrupole moment in Barns
        coefq = 0.75*rayo2/3.1415926
        rfiss = 1.28*aunti - 0.76 + (0.8/aunti)
        rfiss = rfiss*(1.+0.0007*temp*temp)
        rayo3 = rayo2*rayon
        r1 = raya1
        r2 = raya2
        beta = r2/r1
        beta2 = beta*beta
        r1pr2 = r1 + r2
        qupi1 = as*(1.-rkas*z1i2)/(ro*ro)
        qupi2 = as*(1.-rkas*z2i2)/(ro*ro)
        qupig = (qupi1*qupi2)**.5
        z1z2 = z1*z2
        z2sua = z0*z0/a0
        bdif = 0.99*(1.+0.009*temp*temp)
        cosep = .5*bdif*bdif*((r1 + r2)/(r1*r2))
        sep = .5*bdif*bdif*((r1 + r2)/(r1*r2))
        eso = as*(1.-rkas*zii)*adeti
        ess = as*a1det*(1.-rkas*z1i2) + as*a2det*(1.-rkas*z2i2)
        eco = (.864*z0*z0)/rayon
        ecinf = .864*(((z1*z1)/(r1)) + ((z2*z2)/(r2)))
        pas = paray*rayon
        evoinf = avol*((1.-1.80*z1i2)*a1 + (1.-1.80*z2i2)*a2)
        evolsp = avol*(1.-1.80*zii)*a0
        difvol = evoinf - evolsp
        eo = eco + eso + evolsp
        eninf = ess + ecinf + evoinf
        qreac = eo - eninf
        !  qexp becomes qreac of the LDM if qexp=0. in the file fusfis.dat
        IF (qexp == 0.0D0) qexp = qreac
        difes = ess - eso*(r1*r1 + r2*r2)/rayo2
        difec = 0.864*(z1*z1/r1 + z2*z2/r2) + 1.44*z1z2/r1pr2 - (z1 + z2) &
        &        *(z1 + z2)*(.864*(r1**5.+r2**5.) + 1.44*((r1*r2)**3.)/r1pr2) &
        &        /(rayo3*rayo3)
        ! disco=discontinuity at the contact point: two body minus one body energies
        disco = difes + difec + difvol

        z1sa1 = z1/a1
        z2sa2 = z2/a2

        ! TODO
        amo = 0.320*1.2249
        shell = 0.0

        ph1 = 2.*r1*r1
        ph2 = 4.*r1*r1*r1*r1
        ph3 = 2.*ph1
        phs1 = 2.*r2*r2
        phs2 = 4.*r2*r2*r2*r2
        phs3 = 2.*phs1
        denbc = 4.1887902*rayo3*rayo2
        coec = 1.439965*z1*z2
        xuu = r2
        z4 = .5*xuu
        d5 = .4840801*xuu
        c7 = .4180156*xuu
        c8 = .3066857*xuu
        c9 = .1621267*xuu

        do j = nr1, nmax, npas
            if (j <= n1cor) then ! begin one-body shapes (QMS)
                s = (1.*(n1cor - j))/(1.*(n1cor - 1))
                if (j <= 1) s = 0.9995
                if (j < nr1) s = 0.9995
                if (j == n1cor) s = 0.009
                s2 = s*s
                s4 = s2*s2
                fs = s2 + (1.-s2)*beta2
                sprim = s/sqrt(fs)
                sp2 = sprim*sprim
                sp4 = sp2*sp2
                c3 = rayo3*(8./3.)/(v(s) + ((s/sprim)**3.)*v(sprim))
                c = c3**.3333333333
                c2 = c*c
                c2s2 = c2*s2
                cp2 = c2s2/sp2
                c5 = c3*c2
                cp = c*s/sprim
                cp3 = cp2*cp
                cp4 = cp2*cp2
                c4 = c2*c2
                aco = c*s
                ac2 = aco*aco
                ac4 = ac2*ac2
                ac6 = ac4*ac2
                cpac = cp2 - ac2
                c2ac = c2 - ac2
                d2 = .25*ac2*ac2/cpac
                d = sqrt(d2)
                csd = cp/d
                rac = sqrt(cpac)
                arc = log(csd + sqrt(csd*csd + 1.))
                vo1 = 0.5*ac2*cp - cp2*cp/3.+0.5*rac*(cp*sqrt(cp2 + d2) + d2*arc) &
               & - 4.*r2*r2*r2/3.
                xli = 0.
                xri = cp
                iend = 10000
                eps = 1.d-4
                ier = 0

                if (abs(a1 - a2) < 0.01) then
                    rai(j) = c4*(1.0 + s2 + s4)/(4.0*rayo3)
                    rat12 = 0.75*rayon
                    cpt12 = 1.0
                    rcent = rai(j)
                else
                    xl1 = xli
                    xr = xri
                    x = xl1
                    tol = x
                    f = -.5*ac2*x + x*x*x/3.-.5*rac*(x*sqrt(x*x + d2) + d2* &
                    & log(x/d + sqrt(x*x/d2 + 1.))) + vo1
                    if (f == 0) goto 216
                    fl = f
                    x = xr
                    tol = x
                    f = -.5*ac2*x + x*x*x/3.-.5*rac*(x*sqrt(x*x + d2) + d2* &
                    & log(x/d + sqrt(x*x/d2 + 1.))) + vo1
                    if (f == 0) goto 216
                    fr = f
                    if (sign(1.d0, Fl) + sign(1.d0, fr) /= 0) then
                        ier = 2
                        goto 216
                    else
                        i = 0
                        tolf = 100.*eps
                    end if
204                 i = i + 1
                    do m = 1, iend
                        x = .5*(xl1 + xr)
                        tol = x
                        f = -.5*ac2*x + x*x*x/3.-.5*rac*(x*sqrt(x*x + d2) + d2* &
                        & log(x/d + sqrt(x*x/d2 + 1.))) + vo1
                        if (f == 0) goto 216
                        if (sign(1.d0, f) + sign(1.d0, fr) == 0) then
                            tol = xl1
                            xl1 = xr
                            xr = tol
                            tol = fl
                            fl = fr
                            fr = tol
                        end if
                        tol = f - fl
                        za = f*tol
                        za = za + za
                        if (za < fr*(fr - fl)) then
                            if (i <= iend) goto 217
                        end if
                        xr = x
                        fr = f
                        tol = eps
                        za = dabs(xr)
                        if (za > 1.) tol = tol*za
                        if (dabs(xr - xl1) <= tol) then
                            if (dabs(fr - fl) <= tolf) goto 214
                        end if
                    end do
                    ier = 1
214                 if (dabs(fr) > dabs(fl)) then
                        x = xl1
                        f = fl
                    end if
                    go to 216
217                 za = fr - f
                    dx = (x - xl1)*fl*(1.+f*(za - tol)/(za*(fr - fl)))/tol
                    xm = x
                    fm = f
                    x = xl1 - dx
                    tol = x
                    f = -.5*ac2*x + x*x*x/3.-.5*rac*(x*sqrt(x*x + d2) + d2* &
                   & log(x/d + sqrt(x*x/d2 + 1.))) + vo1
                    if (f /= 0) then
                        tol = eps
                        za = dabs(x)
                        if (za > 1.) tol = tol*za
                        if (dabs(dx) <= tol) then
                            if (dabs(f) <= tolf) goto 216
                        end if
                        if (sign(1.d0, f) + sign(1.d0, fl) == 0) then
                            xr = x
                            fr = f
                        else
                            xl1 = x
                            fl = f
                            xr = xm
                            fr = fm
                        end if
                        goto 204
                    end if
216                 zv = x
                    zv2 = zv*zv
                    zv4 = zv2*zv2
                    cpzv = cp2 - zv2
                    czv = c2 - zv2
                    cp4z = cp4 - zv4
                    c4z = c4 - zv4
                    sqzv = (zv2*cpac + .25*ac4)**1.5
                    sqcp = (cp2*cpac + .25*ac4)**1.5
                    sqc = (c2*c2ac + .25*ac4)**1.5
                    zz1 = (.25*ac2*czv - .25*c4z - uti*(sqzv - .125*ac6) &
                    &   /cpac + uti*(sqc - .125*ac6)/c2ac)/(4.*uti*r1*r1*r1)
                    zz2 = (.25*ac2*cpzv - .25*cp4z + uti*(sqcp - sqzv)/cpac) &
                    &   /(4.*uti*r2*r2*r2)
                    rcent = zz1 + zz2
                    if (cpt12 .eq. 0.d0) then
                        rat12 = rcent
                        cpt12 = 1.d0
                    end if
                    r13 = r1*r1*r1
                    r23 = r2*r2*r2
                    xcent = (r13*zz1 - r23*zz2)/(r13 + r23)
                    rai(j) = rcent
                end if

                sq1 = dsqrt(1.-s4)
                sq2 = sq1/s2
                sqp1 = dsqrt(1.-sp4)
                sqp2 = sqp1/sp2
                ssur1 = s4/sq1
                ssup1 = sp4/sqp1
                bs = .25*c2*(1.+ssur1*argsh(sq2) + (s2/sp2) &
                            &  *(1.+ssup1*argsh(sqp2)))/rayo2
                es = eso*bs

                !----------------------------------------------------------
                ! begin the damping of shell and pairing for one-body shapes

                npasin = 40
                pisur2 = pi/2.
                pasin = pisur2/(npasin - 1)

                do intrx = 1, npasin
                    tetax = (intrx - 1)*pasin
                    sinox = sin(tetax)
                    sin2x = sinox*sinox
                    cos2x = 1.0 - sin2x
                    r2x = s*s*sin2x + cos2x
                    wx = sqrt(r2x)*c
                    rrx = s4*sin2x + cos2x
                    rrr = sqrt(rrx)*c*c
                    rintx(intrx) = (wx - rayon)*(wx - rayon)*rrr*sinox
                end do

                call qsf(pasin, rintx, devr, npasin)
                devr1 = devr(npasin)

                do intrx = 1, npasin
                    tetax = (intrx - 1)*pasin + pisur2
                    sinox = sin(tetax)
                    sin2x = sinox*sinox
                    cos2x = 1.0 - sin2x
                    r2x = sp2*sin2x + cos2x
                    wx = sqrt(r2x)*cp
                    rrx = sp4*sin2x + cos2x
                    rrr = sqrt(rrx)*cp*cp
                    rintx(intrx) = (wx - rayon)*(wx - rayon)*rrr*sinox
                end do

                call qsf(pasin, rintx, devr, npasin)
                devr2 = devr(npasin)

                ecart = (0.5*(devr1 + devr2))/(rayon*rayon*es/eso)
                ecart = ecart/(amo*amo)
                amort = (1.0 - 2.0*ecart)*dexp(-ecart)
                couch = shell*amort

                ! end the damping of shell and pairing for one-body shapes
                !----------------------------------------------------------

                !-------------------------------------------------------
                ! begin of the coulomb energy for one-body shapes

                s2m1 = s2 - 1.
                sp2m1 = sp2 - 1.
                do jj = 1, nh1
                    ooo = 0.0001 + (jj - 1)*hv
                    siooo = dsin(ooo)
                    coooo = dcos(ooo)
                    rcc = c2s2*siooo*siooo
                    co2oo = coooo*coooo
                    if (ooo <= 1.570796327) then
                        r2j = rcc + c2*co2oo
                    else
                        r2j = rcc + cp2*co2oo
                    end if
                    rjj = dsqrt(r2j)
                    akr = 4.0*rjj*siooo
                    do k = 1, nh2
                        oi = 0.0003 + (k - 1)*hy
                        sinoi = dsin(oi)
                        cosoi = dcos(oi)
                        cos2 = cosoi*cosoi
                        sin2 = sinoi*sinoi
                        css = c2s2*sin2
                        ccc = c2*cos2
                        cpc = cp2*cos2
                        ooooi = ooo + oi
                        coooi = dcos(ooooi)
                        siooi = dsin(ooooi)
                        cor = rjj*coooi*sinoi
                        cor1 = rjj*sin2*cosoi*siooi
                        de2 = 2.0*rjj*coooi
                        if (oi <= 1.570796327) then
                            r2k = css + ccc
                            rkk = dsqrt(r2k)
                            cor2 = rkk*sinoi
                            coefk = cor*r2k - r2k*cor2 + s2m1*cor1*c2
                            deooi = r2j + r2k - de2*rkk
                            coefd = -0.5*deooi*(-cor2 + (s2m1*ccc*sinoi)/rkk)
                        else
                            r2k = css + cpc
                            rkk = dsqrt(r2k)
                            cor2 = rkk*sinoi
                            deooi = r2j + r2k - de2*rkk
                            coefk = cor*r2k - cor2*r2k + sp2m1*cor1*cp2
                            coefd = -0.5*deooi*(-cor2 + (sp2m1*sinoi*cpc)/rkk)
                        end if
                        ak2 = (akr*rkk*sinoi)/deooi
                        aak = dsqrt(ak2)
                        cefd2 = coefd*ak2
                        y1(k) = (coefk*cei1(aak) + cefd2*cei2(aak))/(dsqrt(deooi)*(-1.0))
                    end do
                    call qsf(hy, y1, vo, nh2)
                    voo(jj) = vo(nh2)*siooo*r2j*rjj
                end do
                call qsf(hv, voo, bcc, nh1)
                bc = bcc(nh1)/denbc
                ec = eco*bc

                ! end of the coulomb energy for one-body shapes
                ! ------------------------------------------------------------

                ! ------------------------------------------------------------
                ! begin of the proximity energy for one-body shapes

                if (s - 0.7071d0 > 0.0d0 .or. sprim - 0.7071d0 > 0.0d0) then
                    en = 0.0d0
                else
                    w1 = 0.5d0*c2
                    w2 = 0.25d0*c2*c2
                    w3 = (1.0d0 - s2)*c2
                    ws1 = w1*s2/sp2
                    ws2 = w2*s4/sp4
                    ws3 = (1.0d0 - sp2)*cp2
                    xl = s*c
                    xu = 0.5d0*xl/(sprim*dsqrt(1.0d0 - sp2))
                    b = xu - xl
                    a = 0.5d0*(xu + xl)
                    e = 0.4840801d0*b
                    y = 0.04063719d0*(phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a + e) &
                    & + phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a - e))
                    e = 0.4180156d0*b
                    y = y + 0.09032408d0*(phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a + e) &
                    & + phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a - e))
                    e = 0.3066857d0*b
                    y = y + 0.1303053d0*(phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a + e) &
                    & + phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a - e))
                    e = 0.1621267d0*b
                    y = y + 0.1561735d0*(phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a + e) &
                    & + phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a - e))
                    y = B*(y + 0.1651197d0*phi(bdif, sep, w1, w2, w3, ws1, ws2, ws3, a))
                    en = qupig*y
                end if

                ! end of the proximity energy for one-body shapes
                ! -------------------------------------------------------------

                ! -------------------------------------------------------------
                ! taking into account the discontinuity and qexp at the contact point
                etota = disco*(rai(j) - rai(1))/(r1pr2 - rai(1))
                etota = etota + (qexp - qreac)*(-rai(j) + r1pr2)/(r1pr2 - rai(1))
                etota = etota + en + ec + es + evolsp - (ess + ecinf + evoinf) - (qexp - qreac)
                ! end taking into account the discontinuity and qexp at the contact point
                ! -------------------------------------------------------------

                ! -------------------------------------------------------------
                ! begin the moments and rotation energy for one-body shapes
                s6 = s2*s4
                rayo5 = rayo2*rayo3
                if (abs(a1 - a2) < 0.01) then
                    xxmoq = argsh(sqrt((1.0 - s2)/s2))/sqrt(1.0 - s2)
                    xmom2 = 112.0/s2 + 8.0 + 30.0*s2 - 135.0*s4 + (120.0*s4 - 135.0*s6)*xxmoq
                    xmom2 = xmom2*c3*c2*s2/(512.0*(1.0 - s2)*rayo5)
                else
                    xmo1 = argsh(2.0*sqrt(1.0 - s2)/s2)*s2/2.0/sqrt(1.0 - s2)
                    xmo1 = xmo1*s4*(8.0 - 9.0*s2)/128.0/(1.0 - s2) - 1.0/15.0 + s4/8.0
                    xmo1 = c3*c2*(xmo1 + (2.0 - s2)*(8.0 - 7.0*s4)/128.0/(1.0 - s2))
                    xmo2 = argsh(2.0*sqrt(1.0 - sp2)/sp2)*sp2/2.0/sqrt(1.0 - sp2)
                    xmo2 = xmo2*sp4*(8.0 - 9.0*sp2)/128.0/(1.0 - sp2) - 1.0/15.0 + sp4/8.0
                    xmo2 = cp2*cp3*(xmo2 + (2.0 - sp2)*(8.0 - 7.0*sp4)/128.0/(1.0 - sp2))
                    xmo1 = xmo1*1.875/rayon**5
                    xmo2 = xmo2*1.875/rayon**5
                    xmom2 = xmo1 + xmo2 - 2.5*xcent**2/rayon**2
                end if
                erot = 51.8325*xlmom*(xlmom + 1.)/(a0*rayo2*xmom2)
                ! 51.8325=0.5*hbar2/m0/0.4 with m0=939.0327 MeV
                ! end the moments and rotation energy
                ! ----------------------------------------------------------------

                ! store the total energy
                etota = etota + erot + couch
                etot(j) = etota
                ! end store the total energy

                if (abs(a1 - a2) < 0.01) then
                    xxmoq = argsh(sqrt((1.0 - s2)/s2))/sqrt(1.0 - s2)
                    quadr = 16.0/s2 - 8.0 - 14.0*s2 + 15.0*s4 - (24.0*s4 - 15.0*s2*s4)*xxmoq
                    quadr = 3.141593*quadr*c3*c2*s2/(96.0*(1.0 - s2)*rayo5)
                    q2 = 0.6*xmom2 - 0.0596831*quadr
                    q2 = q2*rayo2
                else
                    quadr = 1.0
                    q2 = 1.0
                end if

            end if ! end one-body shapes (QMS)

            ! end one-body shapes
            ! -------------------------------------------------------------

            ! -------------------------------------------------------------
            ! begin two-body shapes

            if (j > n1cor) then ! begin calculation for two-body shapes
                s = 0.0
                couch = 0.0
                sprim = 0.0
                c = 2.*r1
                ier = 0
                rcent = r1pr2 + (j - n1cor)*pas
                rai(j) = rcent
                ec = ecinf + coec/rcent
                ! 两种计算离心式方法，哪个对？
                ! method 1:
                ecer = 20.8000*llll*(llll + 1.000)/((a1*a2/a0)*rcent*rcent)

                ! method 2:
                if (abs(a1 - a2) < 0.001) then
                    xmom2 = .62996052 + (.625*rcent*rcent)/(rayo2)
                else
                    r13 = r1*r1*r1
                    r23 = r2*r2*r2
                    xmom2 = 0.4*(r1**5 + r2**5) &
                    & + r13*r23*rcent*rcent/(r13 + r23)
                    xmom2 = xmom2/(0.4*rayo2*(r13 + r23))
                end if
                ecer = xlmom*(xlmom + 1.)*51.8325/(a0*rayo2)/xmom2

                es = ess
                sep = (j - n1cor)*pas + cosep
                u = 0.04063719*(phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 + d5) &
                & + phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 - d5))
                u = u + 0.09032408*(phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 + c7) &
                & + phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 - c7))
                u = u + 0.1303053*(phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 + c8) &
                & + phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 - c8))
                u = u + 0.1561735*(phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 + c9) &
                & + phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4 - c9))
                u = xuu*(u + 0.1651197*phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, z4))
                en = qupig*u

                if (abs(a1 - a2) < 0.01) then
                    xmom2 = 0.62996052 + (0.625*rcent*rcent)/rayo2
                else
                    r13 = r1*r1*r1
                    r23 = r2*r2*r2
                    xmom2 = 0.4*(r1**5 + r2**5) + r13*r23*rcent*rcent/(r13 + r23)
                    xmom2 = xmom2/(0.4*rayo2*(r13 + r23))
                end if

                erot = xlmom*(xlmom + 1.)*51.8325/(a0*rayo2)/xmom2
                etota = ec + en + ecer - ecinf + erot + couch - (qexp - qreac)
                etot(j) = etota

                quadr = 2.094395*rcent*rcent/rayo2
                q2 = 0.37797632*rayo2 + 0.25*rcent*rcent

                ! end two body
            end if ! end calculation for two-body shapes

        end do

        ! The essential characteristics of the barrier are given by
        ! the following lines (the curve has been smoothed before)
        ! The derivatives of the potential are given for dynamics

        do nema = nr1p1, nmax1
            croi = (etot(nema) - etot(nema - 1))*(etot(nema + 1) - etot(nema))
            if (croi <= 0) then
                ebarr = etot(nema) + (qexp - qreac)
                ebafi = ebarr + eninf - eo - (qexp - qreac)
                rbarr = rai(nema)
                wh = (etot(nema + 1) + etot(nema - 1) - 2.0d0*etot(nema))/(pas*pas)
                if (wh <= 0) then
                    wh = (dsqrt(-wh)*6.4421d0)/((a1*a2/a0)**0.5d0)
                    nem = nema
                end if
            end if
        end do

        econt = etot(n1cor + 1) + (qexp - qreac)
        step = 0.1
        jr = 1
        rt(jr) = (int(rai(nm)/step) + 1)*step

        do while (rt(jr) - rai(nmax - 1) <= 0)
            do ir = nm, nmax
                if (rt(jr) - rai(ir) <= 0) exit
            end do
            d9 = dabs(rt(jr) - rai(ir))
            d1 = dabs(rt(jr) - rai(ir - 1))
            if (d9 > d1) ir = ir - 1
            h1 = rai(ir) - rai(ir - 1)
            h2 = rai(ir + 1) - rai(ir)
            h1ph2 = h1 + h2
            fm1 = etot(ir - 1)
            fo = etot(ir)
            fp1 = etot(ir + 1)
            axxo = (rt(jr) - rai(ir))*2.0
            dvj = fm1*(axxo - h2)/h1/h1ph2 &
            &   + fp1*(axxo + h1)/h2/h1ph2 &
            &   + fo*(-axxo + h2 - h1)/h1/h2
            erego(jr) = dvj*(rai(ir) - rt(jr))*(-1.0) + fo
            jr = jr + 1
            rt(jr) = rt(jr - 1) + step
        end do

        kr = jr - 1

        do jr = 1, kr
            if (nr1 - n1cor + 5 > 0) exit
            if (rt(jr) - r1pr2 + 0.1 > 0) then
                jrt = jr
                penta = (erego(jr + 5) - erego(jr + 4))/h
                pentb = (erego(jr - 1) - erego(jr - 2))/h
                di4m1 = erego(jr + 4) - erego(jr - 1)
                dpent = penta - pentb
                paeb = 2.5*(penta + pentb)
                erego(jr) = erego(jr - 1) + (pentb + 0.166666666*dpent)*di4m1/paeb
                erego(jr + 1) = erego(jr) + (pentb + 0.3333333*dpent)*di4m1/paeb
                erego(jr + 2) = erego(jr + 1) + (pentb + 0.5*dpent)*di4m1/paeb
                erego(jr + 3) = erego(jr + 2) + (pentb + 0.6666666*dpent)*di4m1/paeb
                exit
            end if
        end do

        krm1 = kr - 1

        do jd = 2, krm1
            deri(jd) = (erego(jd + 1) - erego(jd - 1))/(2.0*h)
        end do

        ! tabulate the fusion and fission curve
        erefi(1:krm1) = erego(1:krm1) - qreac
        erego(1:krm1) = erego(1:krm1) + (qexp - qreac)

        ! print fusion and fission curve
        if (opt2 /= 0) then
            open (11, file='fusbar.dat', status='replace')
            write (*, 774)
            write (*, 410) (rt(jr), erego(jr), deri(jr), jr=2, krm1)
            write (11, 774)
            write (11, 410) (rt(jr), erego(jr), deri(jr), jr=2, krm1)
            close (11)
        else
            open (12, file='fisbar.dat', status='replace')
            write (*, 874)
            write (*, 410) (rt(jr), erefi(jr), deri(jr), jr=2, krm1)
            write (12, 874)
            write (12, 410) (rt(jr), erefi(jr), deri(jr), jr=2, krm1)
        end if

774     format(/, ' ENERGY RELATIVELY TO THE INFINITY: FUSION BARRIER',/)
874     format(' ENERGY RELATIVELY TO THE GROUND STATE: FISSION BARRIER', /, &
& t10, 'r', t20, 'e', t28, 'der')
410     format( f10.4, f10.4, f10.4)

        ! ----------------------------------------------------------
        ! half-life determination
        do i = 1, krm1 - 1
            if (rt(i) <= r1pr2 .and. rt(i + 1) > r1pr2) then
                ncont = i
            end if
        end do

! 遍历erefi和rt数组，查找潜在的最小值
        if (opt11 == 0 .and. opt3 == 0) then
            redmas = (a1*a2)/a0
            coeff1 = 1.d0/20.9d0
            nminim = 0
            n111 = krm1 - 1
            n222 = krm1 - 2
            if (erefi(1) < erefi(2)) then
                nminim = nminim + 1
                nmin12(nminim) = 1
            end if
            do ij = 2, n111
                if (erefi(ij - 1) > erefi(ij) .and. erefi(ij) < erefi(ij + 1) &
                & .and. rt(ij) < r1pr2) then
                    nminim = nminim + 1
                    nmin12(nminim) = ij
                end if
            end do
            if (nminim > 0) then
                nabsmi = nmin12(1)
                do ji = 1, nminim
                    ikk = nmin12(ji)
                    ! In the following instruction lt corresponds to the tunneling
                    ! from the deepest well and gt would correspond to the highest well
                    ! i.e the isomeric state
                    ! THE ENERGY REFERENCE IS THE SPHERE ENERGY
                    ! AND THE INTEGRATION IS DONE FROM R1+R2     EREF=0.
                    if (opt12 == 0) then
                        if (erefi(ikk) < erefi(nabsmi)) nabsmi = ikk
                    else
                        if (erefi(ikk) > erefi(nabsmi)) nabsmi = ikk
                    end if
                end do
                write (*, 10003) rt(nabsmi), erefi(nabsmi)
10003           format(/, ' MINIMUM CHARACTERISTICS:', ' R=', f14.3, ' E=', f12.3,/)
            else
                write (*, 10001)
10001           format(' THERE IS NO POTENTIAL ENERGY MINIMUM',/)
            end if
        end if

        nturn = 0
        eref = 0.0
        do ii = nabsmi, n222 - 1
            if (erefi(ii) >= eref .and. erefi(ii + 1) <= eref .and. rt(ii) > r1pr2) then
                nturn = ii
                r1t = rt(nturn)
                r2t = rt(nturn + 1)
                vpot1 = erefi(nturn)
                vpot2 = erefi(nturn + 1)
                rout = (eref - vpot1)/(vpot2 - vpot1)*(r2t - r1t) + r1t
                routth = 0.7200*z2*z1/qexp &
                & + sqrt(routth*routth + 20.8*llll*(llll + 1.000)/qexp)
                actn = (sqrt(coeff1*bmas(rat12, rayon, rt(ncont + 1), r1pr2, redmas) &
                & *abs((erefi(ncont + 1) - erefi(ncont))))*(rt(ncont + 1) - rt(ncont)) &
                & + sqrt(coeff1*bmas(rat12, rayon, r1t, r1pr2, redmas) &
                & *abs((vpot1 - eref)))*(rout - r1t))*0.5d0
                do l = ncont + 1, nturn - 1
                    ll = l + 1
                    actn = actn + (sqrt(coeff1*bmas(rat12, rayon, rt(ll), r1pr2, redmas) &
                    & *abs(erefi(ll) - eref)) &
                    & + sqrt(coeff1*bmas(rat12, rayon, rt(l), r1pr2, redmas) &
                    & *abs(erefi(l) - eref)))*(rt(ll) - rt(l))*0.5d0
                end do
                write (*, "(4F12.6, I5)") rt(nturn), rt(nturn + 1), rout, routth, nturn
                exit
            end if
        end do

        if (nturn == 0) stop 'THERE IS NO TURNING POINT!'

        pp12 = exp(-2.0d0*actn)
        freq = 5.0d21
        xlamb = freq*pp12
        t12 = log(2.0d0)/xlamb
        xlog10 = log10(t12)

        ! output
        t12_out = t12

        write (*, "('action = ', f14.6, /, 't1/2 = ', e14.6, ' s', /, &
        & 'log10(t12) = ', f14.6, /, 'pp12 = ', e14.6, /)") &
        & actn, t12, xlog10, pp12

        ! ----------------------------------------------------------------------------
        ! begin fusion cross section with WONG FORMULA

        if (opt10 == 1) then
            tzhw = 0.0
            do jm = 2, krm1
                prodd = (erego(jm) - erego(jm - 1))*(erego(jm + 1) - erego(jm))
                concd = erego(jm) - erego(jm - 1)
                if ((prodd <= 0.0) .and. (concd >= 0.0)) then
                    tzhw = 1.0
                    ! zhw : second derivative of E and after hbaromega
                    zhw = (erego(jm + 1) + erego(jm - 1) - 2.0*erego(jm))/(1.0*HH)
                    zhw = (sqrt(abs(zhw))*6.4421)/((A1*A2/A0)**0.5)
                end if
            end do
            if (tzhw <= 0.0) then
                eee = ebarr - 11.0
                write (*, 1753) zhw
                ! print the cross section for each energy level
                do ijk = 1, 100
                    eee = eee + 1.0
                    auxse = 1.0 + exp(depi*(eee - ebarr)/zhw)
                    secte = (rbarr*rbarr*zhw/(2.0*eee))*log(auxse)*10.0
                    write (*, 5768) eee, secte
                end do
            end if
        end if

1753    format(/, ' cross section (WONG FORMULA valid for medium E):', &
                 & ' HW(MeV)=', f6.4)
5768    format(' E=', f9.5, ' section efficace(mb)=', e19.5)

        ! end fusion cross section with WONG FORMULA
        ! -----------------------------------------------------------------------

        !------------------------------------------------------------------------
        ! exact calculation of the fusion cross section

        if (opt5 == 1) then
            do ie = 1, nsec1, nsec2
                e = 1.*idint(ebarr - 15.+1.*ie)
                write (6, 1600) e, 1/e
1600            format(/, 1x, ' ENERGY(MASS CENTRE):', F6.2, ' 1/E=', F7.4)
                pila2 = 65.188*a0/(a1*a2*e)
                tvar1 = 1.0d0
                do l = 1, 219
                    ll5 = l - 1
                    do lr = 1, kr
                        et1(lr) = erego(lr) + (20.75*l*a0*(l - 1)) &
                        & /(a1*a2*rt(lr)**2)
                    end do
                    test = -1.
                    jrm1 = kr - 1
                    do jm = 2, jrm1
                        prod = (et1(jm) - et1(jm - 1))*(et1(jm + 1) - et1(jm))
                        conca = et1(jm) - et1(jm - 1)
                        if ((prod <= 0.) .and. (conca >= 0.)) then
                            test = 1.
                            lmax = l
                            elmax = et1(jm)
                            hw = (et1(jm + 1) + et1(jm - 1) - 2.*et1(jm))/(1.*hh)
                            hw = (sqrt(dabs(hw))*6.4421)/((a1*a2/a0)**.5)
                            rol = rt(jm)
                            jol = jm
                        end if
                    end do

                    if (test <= 0) exit
                    if (e < elmax) then
                        do lex = 1, jrm1
                            if (rt(lex) > rai(irai)) then
                                chera = (et1(lex) - e)*(et1(lex + 1) - e)
                                if (chera < 0.) then
                                    l1 = lex
                                    ra = rt(lex)
                                    if (ra <= rol) exit
                                    lmax = lmax - 1
                                    goto 1619
                                end if
                            end if
                        end do

                        do lex = l1 + 1, jrm1
                            cherb = (et1(lex) - e)*(et1(lex + 1) - e)
                            if (cherb <= 0.) then
                                rb = rt(lex)
                                l2 = lex
                                exit
                            end if
                        end do
                        ndim = l2 - l1 + 1
                        do lrr = 1, ndim
                            lvl = lrr + l1 - 1
                            evme = abs(et1(lvl) - e)
                            et1(lrr) = sqrt(evme)
                        end do
                        call qsf(h, et1, zz, ndim)
                        rkl = .219527*((a1*a2/a0)**.5)*zz(ndim)
                        tl(l) = 1./(1.+dexp(2.*rkl))
                        if (tl(l) <= 1.0d-8) exit
                    else
                        deshw = depi/hw
                        tl(l) = 1./(1.+dexp(deshw*(elmax - e)))
                        rb = 0.
                        ra = 0.
                        ndim = 0
                        if (tl(l) <= 1.0d-8) exit
                    end if
                    tvar2 = tl(l)
                    vartl = tvar2 - tvar1
                    tvar1 = tl(l)
                    if (vartl > 0.) exit
                    if (opt6 == 0) then
                        write (*, 1695) ll5, elmax, rol, jol, ra, rb, ndim, tl(l), hw
1695                    format(' L=', I3, ' ELMAX=', F8.3, ' ROL=', F6.3, ' JOL=', &
                                                         & I3, ' RA=', F6.3, /, ' RB=', F6.3, ' NDIM=', I3, &
                                                         & ' TL=', F9.7, ' HW=', E11.3)
                    end if

                    l7 = ll5
                    emaxl = elmax
                    rrol = rol
                    jjol = jol
                    rra = ra
                    rrb = rb
                    nndim = ndim
                    hww = hw
                    ttl = tl(l)
                end do

                WRITE (*, 1720) l7, emaxl, rrol, jjol, rra, rrb, nndim, ttl, hww
1720            format(' L=', I3, ' ELMAX=', F8.3, ' ROL=', F6.3, ' JOL=', I3, &
                                                    & ' RA=', F6.3, /, '  RB=', F6.3, ' NDIM=', I3, ' TL=', F9.7, &
                                                    & ' HW=', E11.3)
1619            sect = 0.
                do le = 1, lmax
                    sect = sect + pila2*(2.*(le - 1) + 1.)*tl(le)
                end do
                sect = sect*10.
                write (*, 1603) sect, lmax
1603            format(' fusion cross section:', E12.5, ' LCR=', I3)
            end do
        end if

        ! end exact calculation of the fusion cross section
        !------------------------------------------------------

        write (*, *) "here is ok !"

    end subroutine gldm

    !======================================================================
    ! sub and func below
    ! =====================================================================
    function v(x)
        real(8) :: x, v
        v = 1./3.+x*x*0.5 + &
        &   0.25*x*x*x*x*log((2.-x*x + 2.*sqrt(1.-x*x))/x/x)/sqrt(1.-x*x)
    end function v

    function argsh(x)
        real(8) :: x, argsh
        argsh = log(x + sqrt(x*x + 1.))
    end function argsh

    subroutine qsf(h, y, z, ndim)
        implicit none
        integer, intent(in) :: ndim
        real(8), intent(in) :: h
        real(8), dimension(100), intent(in) :: y
        real(8), dimension(100), intent(out) :: z
        integer :: i, l1, l2, l3, l4, l5, l6
        real(8) :: ht, sum1, sum2, aux1, aux2, aux

        ht = 0.33333333*h
        l1 = 1
        l2 = 2
        l3 = 3
        l4 = 4
        l5 = 5
        l6 = 6
        if (ndim > 5) then
            ! preparations of integration loop
            sum1 = y(l2) + y(l2)
            sum1 = sum1 + sum1
            sum1 = ht*(y(l1) + sum1 + y(l3))
            aux1 = y(l4) + y(l4)
            aux1 = aux1 + aux1
            aux1 = sum1 + ht*(y(l3) + aux1 + y(l5))
            aux2 = ht*(y(l1) + 3.875*(y(l2) + y(l5)) + 2.625*(y(l3) + y(l4)) + y(l6))
            sum2 = y(l5) + y(l5)
            sum2 = sum2 + sum2
            sum2 = aux2 - ht*(y(l4) + sum2 + y(l6))
            z(l1) = 0.
            aux = y(l3) + y(l3)
            aux = aux + aux
            z(l2) = sum2 - ht*(y(l2) + aux + y(l4))
            z(l3) = sum1
            z(l4) = sum2
            do i = 7, ndim, 2
                ! integration loop
                sum1 = aux1
                sum2 = aux2
                aux1 = y(i - 1) + y(i - 1)
                aux1 = aux1 + aux1
                aux1 = sum1 + ht*(y(i - 2) + aux1 + y(i))
                z(i - 2) = sum1
                if (i == ndim) then
                    z(ndim - 1) = sum2
                    z(ndim) = aux1
                    return
                end if
                aux2 = y(i) + y(i)
                aux2 = aux2 + aux2
                aux2 = sum2 + ht*(y(i - 1) + aux2 + y(i + 1))
                z(i - 1) = sum2
            end do
            z(ndim - 1) = aux2
            z(ndim) = aux1
        else if (ndim == 4 .or. ndim == 5) then
            ! ndim is equal to 4 or 5
            sum2 = 1.125*ht*(y(l1) + 3.*y(l2) + 3.*y(l3) + y(l4))
            sum1 = y(l2) + y(l2)
            sum1 = sum1 + sum1
            sum1 = ht*(y(l1) + sum1 + y(l3))
            z(l1) = 0.
            aux1 = y(l3) + y(l3)
            aux1 = aux1 + aux1
            z(l2) = sum2 - ht*(y(l2) + aux1 + y(l4))
            if (ndim == 5) then
                aux1 = y(l4) + y(l4)
                aux1 = aux1 + aux1
                z(l5) = sum1 + ht*(y(l3) + aux1 + y(l5))
            end if
            z(l3) = sum1
            z(l4) = sum2
        else if (ndim == 3) then
            ! ndim is equal to 3
            sum1 = ht*(1.25*y(l1) + y(l2) + y(l2) - 0.25*y(l3))
            sum2 = y(l2) + y(l2)
            sum2 = sum2 + sum2
            z(l3) = ht*(y(l1) + sum2 + y(l3))
            z(l1) = 0.
            z(l2) = sum1
        end if
    end subroutine qsf

    ! elliptic integrals to determine the Coulomb energy
    function cei1(am)
        implicit none
        real(8) :: cei1, am, geo, ari, aari, test
        if (1.0d0 <= am*am) then
            cei1 = 1.0d34
            return
        else
            geo = dsqrt(1.0d0 - am*am)
            ari = 1.0d0
            do
                aari = ari
                test = aari*1.0d-4
                ari = ari + geo
                if ((aari - geo) <= test) then
                    cei1 = 3.141592654/ari
                    return
                    exit
                else
                    geo = dsqrt(aari*geo)
                    ari = 0.5d0*ari
                end if
            end do
        end if
    end function cei1

    function cei2(al)
        implicit none
        real(8) :: cei2, al, geom, arit, aa, an, w, aarit
        if (1.0d0 <= al*al) then
            cei2 = 1.0d34
            return
        else
            geom = dsqrt(1.0d0 - al*al)
            arit = 1.0d0
            aa = 0.0d0
            an = 1.0d0
            w = 1.0d0
            do
                w = w + aa*geom
                w = w + w
                aa = an
                aarit = arit
                arit = geom + arit
                an = w/arit + an
                if (aarit - geom <= 1.0d-4*aarit) then
                    cei2 = 0.7853981634d0*an/arit
                    return
                    exit
                else
                    geom = dsqrt(geom*aarit)
                    geom = geom + geom
                end if
            end do
        end if
    end function cei2

    function phi(bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, h) result(y)
        implicit none
        real(8), intent(in) :: bdif, sep, ph1, ph2, ph3, phs1, phs2, phs3, h
        real(8) :: y, h2, aa1, aa2, sqr1, sqr2, dist, dist1, dist2, eta, eta2, eta3, eta4

        h2 = h*h
        aa1 = ph2 - h2*ph3
        aa2 = phs2 - h2*phs3
        sqr1 = dsqrt(aa1)
        sqr2 = dsqrt(aa2)
        dist = (dsqrt(ph1 - h2 - sqr1) + dsqrt(phs1 - h2 - sqr2) + sep)/bdif
        dist2 = dist*dist
        eta = 2.74 - dist
        eta2 = eta*eta
        eta3 = eta2*eta
        eta4 = eta2*eta2

        if (dist <= 1.2311) then
            y = (-1.0d0 + 0.1889d0*dist2)*h
        else
            dist1 = dist - 1.2311
            if (dist1 <= 1.5089) then
                y = (-0.135d0 - 0.1881d0*eta - 0.1581d0*eta2 - 0.01202d0*eta3 + 0.02055d0*eta4)*h
            else
                y = (-6.145d0*dexp(-dist/0.7176d0))*h
            end if
        end if

    end function phi

    function bmas(rat12, rayon, rx, r1pr2, redmas)
        implicit none
        real(8), intent(in) :: rat12, rayon, rx, r1pr2, redmas
        real(8) :: aux, bmas

        if (rx < r1pr2) then
            aux = ((r1pr2 - rx)/(r1pr2 - rat12))**2
            bmas = (16.0d0*17.0d0*aux/15.0d0 &
                   & *dexp(-128.0d0/51.0d0*((rx - rat12)/rayon)) + 1.0d0) &
                  & *redmas
        else
            bmas = redmas
        end if
    end function bmas

end module hello_gldm
