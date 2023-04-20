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
    & aunti, adeti, amo
    real(8) :: beta, beta1, beta2, beta3, bdif
    real(8) :: coefq, cosep
    real(8) :: depi, difvol, difes, difec, disco
    real(8) :: eso, ess, eco, ecinf, evoinf, evolsp, eo, eninf
    real(8) :: h, hh, hv, hy
    real(8) :: pas
    real(8) :: qexp, qvalu, qupi1, qupi2, qupig, qreac
    real(8) :: r1, r2, ro, rkas, raya1, raya2, rfiss, rayon, rayo2, rayo3, &
    & r1pr2
    real(8) :: sep
    real(8) :: t12, tcrit
    real(8) :: uti
    real(8) :: z0, z1, z2, z1i, z1i2, z2i, z2i2, zi, zii, z1z2, z2sua,z1sa1, z2sa2
    integer :: nsec1, nsec2, nh1, nh2, nr1p1, nm, nmax1

    a1 = frag_a1
    z1 = frag_z1
    a2 = frag_a2
    z2 = frag_z2
    qexp = qexp_in
    write(*,*) 'a1, z1, a2, z2, qexp:', a1, z1, a2, z2, qexp


    ! 计算t12
    t12 = 123.
    t12_out = t12
    write(*,*) 't12_out:', t12_out


    ! old main gldm
    a0=a1+a2
    z0=z1+z2

    tcrit = 17.
    qvalu = 0.0

    ! NSEC -7 to -7+NSEC1 STEP NSEC2 to have energy for fusiopn cross sections
    nsec1 = 160
    nsec2 = 8
    nh1 = 25
    nh2 = 29
    hv=3.141470/(1.*(nh1-1)*1.)
    hy=3.141290/(1.*(nh2-1)*1.)
    as = 17.9439*(1.+1.5*temp/tcrit)*((1.-temp/tcrit)**1.5)
    avol =-15.4941*(1.+0.00337*temp*temp)
    ro=1.2249
    rkas=2.6
    uti=1./3.
    nr1p1 = nr1+1
    nm = nr1+1
    nmax1 = nmax-1
    depi = 6.2831853
    h = 0.1
    hh = h*h
    alpha=(a1-a2)/(a1+a2)
    a1unt=a1**.3333333333
    a1det=a1**.6666666666
    a2unt=a2**.3333333333
    a2det=a2**.6666666666
    aunti=a0**.3333333333
    adeti=a0**.6666666666
    z1i=(a1-2.*z1)/a1
    z1i2=z1i*z1i
    z2i=(a2-2.*z2)/a2
    z2i2=z2i*z2i
    zi=(a0-2.*z0)/a0
    zii=zi*zi
    raya1=1.28*a1unt-0.76+(0.8/a1unt)
    raya2=1.28*a2unt-0.76+(0.8/a2unt)
    if (opt4 == 0) then
      beta1 = raya2 / raya1
      beta3 = beta1 * beta1 * beta1
      rayon = 1.28 * aunti - 0.76 + (0.8 / aunti)
      rayon = rayon * (1.0 + 0.0007 * temp * temp)
      raya1 = rayon / (1.0 + beta3) ** 0.33333333
      raya2 = beta1 * raya1
    else
      rayon = (raya1 ** 3.0 + raya2 ** 3.0) ** uti
      rayon = rayon * (1.0 + 0.0007 * temp * temp)
      raya1 = raya1 * (1.0 + 0.0007 * temp * temp)
      raya2 = raya2 * (1.0 + 0.0007 * temp * temp)
    endif
    rayo2 = rayon * rayon

    ! Coefficient to express the quadrupole moment in Barns
    coefq=0.75*rayo2/3.1415926
    rfiss=1.28*aunti-0.76+(0.8/aunti)
    rfiss=rfiss*(1.+0.0007*temp*temp)
    rayo3=rayo2*rayon
    r1=raya1
    r2=raya2
    beta = r2/r1
    beta2=beta*beta
    r1pr2=r1+r2
    qupi1=as*(1.-rkas*z1i2)/(ro*ro)
    qupi2=as*(1.-rkas*z2i2)/(ro*ro)
    qupig=(qupi1*qupi2)**.5
    z1z2=z1*z2
    z2sua=z0*z0/a0
    bdif=0.99*(1.+0.009*temp*temp)
    cosep=.5*bdif*bdif*((r1+r2)/(r1*r2))
    sep=.5*bdif*bdif*((r1+r2)/(r1*r2))
    eso=as*(1.-rkas*zii)*adeti
    ess=as*a1det*(1.-rkas*z1i2)+as*a2det*(1.-rkas*z2i2)
    eco=(.864*z0*z0)/rayon
    ecinf=.864*(((z1*z1)/(r1))+((z2*z2)/(r2)))
    pas=paray*rayon
    evoinf=avol*((1.-1.80*z1i2)*a1+(1.-1.80*z2i2)*a2)
    evolsp=avol*(1.-1.80*zii)*a0
    difvol=evoinf - evolsp
    eo=eco+eso+evolsp
    eninf=ess+ecinf+evoinf
    qreac=eo-eninf
   !  qexp becomes qreac of the LDM if qexp=0. in the file fusfis.dat
    IF ( qexp==0.0D0 ) qexp = qreac
    difes = ess - eso*(r1*r1+r2*r2)/rayo2
    difec = 0.864*(z1*z1/r1+z2*z2/r2) + 1.44*z1z2/r1pr2 - (z1+z2) &
    &        *(z1+z2)*(.864*(r1**5.+r2**5.)+1.44*((r1*r2)**3.)/r1pr2) &
    &        /(rayo3*rayo3)
     disco=difes+difec+difvol             
     z1sa1=z1/a1
     z2sa2=z2/a2


     ! TODO 
     amo = 0.320*1.2249
  


    write(*,*) "here is ok !"


    end subroutine gldm

end module hello_gldm
