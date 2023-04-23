module hello_gldm
    implicit none
    private

    public :: say_hello, gldm

    real(8) :: temp = 0., xlmom = 0., paray = 0.02
    integer :: nr1 = 1, nmax = 500, npas = 1, n1cor = 50
    integer :: opt1 = 1, opt4 = 1, opt5 = 1, opt6 = 1, &
    & opt7 = 1

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
        & aunti, adeti, amo, aco, ac2, ac4, ac6, arc
        real(8) :: beta, beta1, beta2, beta3, bdif
        real(8) :: coefq, cosep, coec, c, c2, c3, c4, c5, c7, c8, c9, &
        & c2s2, cp2, cp, cp3, cp4, cpac, c2ac, csd, cpt12, cpzv, czv, &
        & cp4z, c4z
        real(8) :: depi, difvol, difes, difec, disco, denbc, d5, d2, d, &
        & dx
        real(8) :: eso, ess, eco, ecinf, evoinf, evolsp, eo, eninf, eps
        real(8) :: fs, f, fl, fr, fm
        real(8) :: h, hh, hv, hy
        real(8) :: pas, ph1, ph2, ph3, phs1, phs2, phs3
        real(8) :: qexp, qvalu, qupi1, qupi2, qupig, qreac
        real(8) :: r1, r2, ro, rkas, raya1, raya2, rfiss, rayon, rayo2, rayo3, &
        & r1pr2, rac, rat12, rcent, r13, r23
        real(8) :: sep, s, s2, s4, sprim, sp2, sp4, sqzv, sqcp, sqc
        real(8) :: t12, tcrit, tol, tolf
        real(8) :: uti
        real(8) :: vo1
        real(8) :: x, xuu, xli, xri, xl1, xr, xm, xcent
        real(8) :: z0, z1, z2, z1i, z1i2, z2i, z2i2, zi, zii, z1z2, z2sua, z1sa1, &
        & z2sa2, z4, za, zv, zv2, zv4, zz1, zz2
        real(8), dimension(999) :: rai
        integer :: nsec1, nsec2, nh1, nh2, nr1p1, nm, nmax1, i, j, iend, ier, m

        a2 = frag_a2
        z2 = frag_z2
        z1 = frag_z1
        a1 = frag_a1
        qexp = qexp_in
        write (*, *) 'a1, z1, a2, z2, qexp:', a1, z1, a2, z2, qexp

        ! 计算t12
        t12 = 123.
        t12_out = t12
        write (*, *) 't12_out:', t12_out

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
        disco = difes + difec + difvol
        z1sa1 = z1/a1
        z2sa2 = z2/a2

        ! TODO
        amo = 0.320*1.2249
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
            if (j <= n1cor) then
                ! begin Q M S
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

                ! end Q M S
            end if

            if (j > n1cor) then
                ! begin two body
                s = 0.0
                ! end two body
            end if

        end do

        write (*, *) "here is ok !"

    end subroutine gldm

    function v(x)
        real(8) :: x, v
        v = 1./3.+x*x*0.5 + &
        &   0.25*x*x*x*x*log((2.-x*x + 2.*sqrt(1.-x*x))/x/x)/sqrt(1.-x*x)
    end function v

    function argsh(x)
        real(8) :: x, argsh
        argsh = log(x + sqrt(x*x + 1.))
    end function argsh

end module hello_gldm
