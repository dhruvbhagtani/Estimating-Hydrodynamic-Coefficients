MODULE procedures
IMPLICIT NONE

CONTAINS

FUNCTION inv(A) result(Ainv)
IMPLICIT NONE
!REAL, DIMENSION(:,:), INTENT(IN) :: A
DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN) :: A
DOUBLE PRECISION, DIMENSION(size(A,1),size(A,2)) :: Ainv

DOUBLE PRECISION, DIMENSION(size(A,1)) :: work  ! work array for LAPACK
INTEGER, DIMENSION(size(A,1)) :: ipiv   ! pivot indices
INTEGER :: n, info

! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI

! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(n, n, Ainv, n, ipiv, info)

IF (info /= 0) THEN
STOP 'Matrix is numerically singular!'
END IF

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
CALL DGETRI(n, Ainv, n, ipiv, work, n, info)!

IF (info /= 0) THEN
STOP 'Matrix inversion failed!'
END IF
END FUNCTION inv
END MODULE procedures

PROGRAM TENTH
USE procedures
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793238462
DOUBLE PRECISION, PARAMETER :: g = 9.80665
DOUBLE PRECISION, PARAMETER :: c = 300
DOUBLE PRECISION, PARAMETER :: sts = 45000000
DOUBLE PRECISION, PARAMETER :: m = 1017876
DOUBLE PRECISION, PARAMETER :: R = 2
INTEGER, PARAMETER :: Nt = 3000
DOUBLE PRECISION, PARAMETER :: dt = 0.01
DOUBLE PRECISION, PARAMETER :: dz = 0.01
INTEGER, PARAMETER :: Nz = 6001
DOUBLE PRECISION, PARAMETER :: H = 6
DOUBLE PRECISION, PARAMETER :: rho = 1025
DOUBLE PRECISION, PARAMETER :: D = 1.5
DOUBLE PRECISION, PARAMETER :: TP = 12
DOUBLE PRECISION, PARAMETER :: k = 0.0296
DOUBLE PRECISION, PARAMETER :: dp = 60
DOUBLE PRECISION, PARAMETER :: Num = 1
INTEGER, PARAMETER :: Numi = 1
DOUBLE PRECISION, DIMENSION(1,Numi) :: Cdens, Cmens
DOUBLE PRECISION, DIMENSION(Numi,Nt) :: Cdmat, Cmmat
!DOUBLE PRECISION, DIMENSION(1,Num) :: Cdens, Cmens
DOUBLE PRECISION, PARAMETER :: Cd = 1
DOUBLE PRECISION, PARAMETER :: Cm = 2
DOUBLE PRECISION, PARAMETER :: Cm0 = 5
DOUBLE PRECISION, PARAMETER :: Cd0 = 2.5
DOUBLE PRECISION, DIMENSION(4,4) :: I, Q, F, P, TF, RTM, RTB, TRTB
DOUBLE PRECISION, DIMENSION(4,Nt) :: xt, x
DOUBLE PRECISION, DIMENSION(1,Nt) :: Ic, Dc, Icpg, Dcpg, z, Fd, Fi,vi,vd,value,Cmfinal,Cdfinal
DOUBLE PRECISION, DIMENSION(2,Nt) :: zi
DOUBLE PRECISION, DIMENSION(1,1) :: W
DOUBLE PRECISION :: sumi, sumd, omega, mean, theta, stdev, const, c1m, stsm, temp1, temp2, mo, c1, c2,temp(2),ran,ui
DOUBLE PRECISION, DIMENSION(2,4) :: Hi
DOUBLE PRECISION, DIMENSION(4,2) :: THi, KF1
INTEGER :: tvar, zvar, a, b, j, dnum, enum, var, ks
DOUBLE PRECISION :: xb1,xb2,xb3,xb4,xb5,xb6,xb7,xb8,xb9,us,sumCd,sumCm
DOUBLE PRECISION, DIMENSION(8,8) :: Amat, summat, val, Iden
DOUBLE PRECISION, DIMENSION(2,1) :: Gs, t1s, t2s,temp3
DOUBLE PRECISION, DIMENSION(2,2) :: Fs, Is, Ps, tex, Ri, lol
DOUBLE PRECISION, DIMENSION(2, Nt) :: xs, xts
DOUBLE PRECISION, DIMENSION(1,2) :: Hs
DOUBLE PRECISION, DIMENSION(4,1) :: Gi, t1, t2
DOUBLE PRECISION, DIMENSION(1,4) :: TGi

omega = (2*PI)/TP

a = 1
b = 1
ui = 1
do while(a<=4)
    b = 1
        do while(b<=4)
            if(a == b) then
                I(a,b) = 1
            else
                I(a,b) = 0
            end if
        b = b + 1
    end do
    a = a + 1
end do

a = 1
b = 1
do while  (a<=8)
    b = 1
    do while (b<=8)
        if(a==b) then
            Iden(a,b) = 1
        else
            Iden(a,b) = 0
        end if
        b = b+1
    end do
    a = a+1
end do

Ri(1,1) = 0.7
Ri(1,2) = 0.7
Ri(2,1) = 0.5
Ri(2,2) = 0.9

P(1,1) = 1.1
P(1,2) = 0.5
P(1,3) = 0.3
P(1,4) = 1.9
P(2,1) = 1.3
P(2,2) = 1.7
P(2,3) = 0.3
P(2,4) = 0.7
P(3,1) = 1.5
P(3,2) = 0.2
P(3,3) = 1.4
P(3,4) = 0.5
P(4,1) = 0.3
P(4,2) = 0.2
P(4,3) = 0.7
P(4,4) = 1.2

xt(1,1) = 0.5
xt(2,1) = 0.3
xt(3,1) = Cm
xt(4,1) = Cd

x(1,1) = 0.1
x(2,1) = 0.2
x(3,1) = Cm0
x(4,1) = Cd0

b = 1
mean = 0
stdev = 1
do while(b<=Numi)
    if(stdev < 0.0d0) then
        print*, "Standard Deviation must be +ve"
    else
        call RANDOM_NUMBER(temp)
        ran=(-2.0d0*log(temp(1)))**0.5
        theta = 2.0d0*PI*temp(2)
        const = mean+stdev*ran*sin(theta)*0.1
    end if
    if (abs(const)>0.1*Cd0) then
        const = 0.1 * Cd0
    end if
    Cdens(1,b) = Cd0 + const/2
    Cmens(1,b) = Cm0 + const
    print*,const
    b = b + 1
end do

b = 1
do while(b<=Numi)
    temp2 = 0
    tvar = 1
    do while(tvar<=Nt)
        mo = cos(omega*temp2)
        zvar = 1
        temp1 = 0
        sumi = 0
        sumd = 0
        temp1 = 0
        do while(zvar<=Nz)
            sumi = sumi + ((PI**3*D*D*H*Cmens(1,b)*rho)/(2*TP*TP*sinh(k*dp))) * cosh(k*(dp + temp1)) * sin(-omega*temp2) * dz
            sumd = sumd+((rho*Cdens(1,b)*D*PI*PI*H*H)/(TP*TP*2*sinh(k*dp)*sinh(k*dp)))*(cosh(k*(dp+temp1)))**2 *abs(mo)*mo*dz
            zvar = zvar + 1
            temp1 = temp1 - dz
        end do
        Ic(1,tvar) = sumi
        Dc(1,tvar) = sumd
        z(1,tvar) = Ic(1,tvar) + Dc(1,tvar)
        temp2 = temp2 + dt
        tvar = tvar + 1
    end do
b = b + 1
end do

Is(1,1) = 1
Is(1,2) = 0
Is(2,1) = 0
Is(2,2) = 1
Ps(1,1) = 1.3
Ps(1,2) = 0.3
Ps(2,1) = 0.2
Ps(2,2) = 1.7
xts(1,1) = 0.5
xts(2,1) = 0.3
xs(1,1) = 0.7
xs(2,1) = 0.5

c1m = c/m
stsm = sts/m
c1 = 1 + (c1m*dt)
c2 = 1 + (stsm*dt*dt)/(2)

xb1 = 1 + ((stsm*dt)/(c1m)) - ((stsm)/(c1m**2)) * log(1+c1m*dt)
xb2 = 1
xb3 = ((1)/(c1m)) * log(1+c1m*dt)
xb4 = ((dt/c1m)) - (log(1+c1m*dt))/(c1m**2)
xb5 = (-stsm*dt)/(xb1*c1)
xb6 = 1/c1 - (stsm*dt*xb3)/(c1*xb1)
xb7 = dt/c1 - (stsm*dt * xb4)/(c1*xb1)
Fs(1,1) = 1/xb1
Fs(1,2) = xb3/xb1
Fs(2,1) = xb5
Fs(2,2) = xb6

Gs(1,1) = ((xb4)/(xb1)) * ((z(1,1))/m)
Gs(2,1) = (xb7) * (z(1,1))/m

a = 2
us = 1
do while(a<=Nt)
    t1s(1,1) = xts(1,a-1)
    t1s(2,1) = xts(2,a-1)
    Gs(1,1) = (z(1,a) * xb4)/(xb1*m)
    Gs(2,1) = (z(1,a) * xb7)/m
    t2s = matmul(Fs,t1s) + Gs*us
    xts(1,a) = t2s(1,1)
    xts(2,a) = t2s(2,1)
    a = a+1
end do

F(1,1) = 1/xb1
F(1,2) = xb3/xb1
F(1,3) = ((xb4)/(xb1)) * ((Ic(1,1)))/(m*Cm)
F(1,4) = ((xb4)/(xb1)) * ((Dc(1,1)))/(m*Cd)
F(2,1) = xb5
F(2,2) = xb6
F(2,3) = (xb7) * (Ic(1,1))/(m*Cm)
F(2,4) = (xb7) * (Dc(1,1))/(m*Cd)
F(3,1) = 0
F(3,2) = 0
F(3,3) = 1
F(3,4) = 0
F(4,1) = 0
F(4,2) = 0
F(4,3) = 0
F(4,4) = 1

Gi(1,1) = 1
Gi(2,1) = 0.5
Gi(3,1) = 0.1
Gi(4,1) = 0.2

TGi(1,1) = 1
TGi(1,2) = 0.5
TGi(1,3) = 0.1
TGi(1,4) = 0.2

b = 1
do while(b<=Num)

    x(3,1) = Cmens(1,b)
    x(4,1) = Cdens(1,b)

    temp2 = 0
    tvar = 1
    do while(tvar<=Nt)
        mo = cos(omega*temp2)
        zvar = 1
        temp1 = 0
        sumi = 0
        sumd = 0
        temp1 = 0
        do while(zvar<=Nz)
            sumi = sumi + ((PI**3*D*D*H*Cmens(1,b)*rho)/(2*TP*TP*sinh(k*dp))) * cosh(k*(dp + temp1)) * sin(-omega*temp2) * dz
            sumd = sumd+((rho*Cdens(1,b)*D*PI*PI*H*H)/(TP*TP*2*sinh(k*dp)*sinh(k*dp)))*(cosh(k*(dp+temp1)))**2 *abs(mo)*mo*dz
            zvar = zvar + 1
            temp1 = temp1 - dz
        end do
        Ic(1,tvar) = sumi
        Dc(1,tvar) = sumd
        z(1,tvar) = Ic(1,tvar) + Dc(1,tvar)
        temp2 = temp2 + dt
        tvar = tvar + 1
    end do

    a = 2
    do while (a<=Nt)

        t1(1,1) = xt(1,a-1)
        t1(2,1) = xt(2,a-1)
        t1(3,1) = xt(3,a-1)
        t1(4,1) = xt(4,a-1)

        F(1,3) = ((xb4)/(xb1)) * ((Ic(1,a-1)))/(m*Cm)
        F(1,4) = ((xb4)/(xb1)) * ((Dc(1,a-1)))/(m*Cd)
        F(2,3) = (xb7) * (Ic(1,a-1))/(m*Cm)
        F(2,4) = (xb7) * (Dc(1,a-1))/(m*Cd)

        t2 = matmul(F,t1) !+ Gi*ui
        xt(1,a) = t2(1,1)
        xt(2,a) = t2(2,1)
        xt(3,a) = t2(3,1)
        xt(4,a) = t2(4,1)
        a = a+1
    end do

!The forces we have calculated are true values, i.e.,they don't have any noise induced. However, the field results will always have error.

    mean = 0
    stdev = 1
    a = 1
    do while (a<=Nt)
        if(stdev < 0.0d0) then
            print*, "Standard Deviation must be +ve"
        else
            call RANDOM_NUMBER(temp)
            ran=(-2.0d0*log(temp(1)))**0.5
            theta = 2.0d0*PI*temp(2)
            const = mean+stdev*ran*sin(theta)
        end if
        vi(1,a) = (const) * 1000
        if(abs(vi(1,a))> abs(0.3*Ic(1,a))) then
            vi(1,a) = 0.3 * Ic(1,a)
        end if
        vd(1,a) = (const) * 1000
        if(abs(vi(1,a))> abs(0.3*Dc(1,a))) then
            vd(1,a) = 0.3 * Dc(1,a)
        end if
        zi(2,a) = xts(1,a) + (vd(1,a))/100000
        a = a + 1
    end do
    Icpg = Ic + vi
    Dcpg = Dc + vd

    Hi(1,1) = 0
    Hi(1,2) = 0
    Hi(1,3) = Ic(1,1)/2
    Hi(1,4) = Dc(1,1)
    Hi(2,1) = 1
    Hi(2,2) = 0
    Hi(2,3) = 0
    Hi(2,4) = 0

    THi(1,1) = Hi(1,1)
    THi(2,1) = Hi(1,2)
    THi(3,1) = Hi(1,3)
    THi(4,1) = Hi(1,4)
    THi(1,2) = Hi(2,1)
    THi(2,2) = Hi(2,2)
    THi(3,2) = Hi(2,3)
    THi(4,2) = Hi(2,4)

    a = 1
    do while(a<=Nt)
        zi(1,a) = Dcpg(1,a) + Icpg(1,a)
        a = a+1
    end do

    !W = 10000
    a = 2
    do while (a<=Nt)
        if(stdev < 0.0d0) then
            print*, "Standard Deviation must be +ve"
        else
            call RANDOM_NUMBER(temp)
            ran=(-2.0d0*log(temp(1)))**0.5
            theta = 2.0d0*PI*temp(2)
            const = mean+stdev*ran*sin(theta)
        end if
        if(const>1) then
            const = 1
        end if
        ui = sqrt(R) * (const) * 100
        W = 10000 * sin(0.5 * a)
        F(1,3) = ((xb4)/(xb1)) * ((Ic(1,a-1)))/(m*Cmens(1,b))
        F(1,4) = ((xb4)/(xb1)) * ((Dc(1,a-1)))/(m*Cdens(1,b))
        F(2,3) = (xb7) * (Ic(1,a-1))/(m*Cm)
        F(2,4) = (xb7) * (Dc(1,a-1))/(m*Cd)

        TF(1,1) = F(1,1)
        TF(1,2) = F(2,1)
        TF(1,3) = F(3,1)
        TF(1,4) = F(4,1)
        TF(2,1) = F(1,2)
        TF(2,2) = F(2,2)
        TF(2,3) = F(3,2)
        TF(2,4) = F(4,2)
        TF(3,1) = F(1,3)
        TF(3,2) = F(2,3)
        TF(3,3) = F(3,3)
        TF(3,4) = F(4,3)
        TF(4,1) = F(1,4)
        TF(4,2) = F(2,4)
        TF(4,3) = F(3,4)
        TF(4,4) = F(4,4)

        t1(1,1) = x(1, a-1)
        t1(2,1) = x(2, a-1)
        t1(3,1) = x(3, a-1)
        t1(4,1) = x(4, a-1)
        t2 = matmul(F,t1) + Gi*ui
        x(1,a) = t2(1,1)
        x(2,a) = t2(2,1)
        x(3,a) = t2(3,1)
        x(4,a) = t2(4,1)

        Hi(1,1) = 0
        Hi(1,2) = 0
        Hi(1,3) = Icpg(1,a)/2
        Hi(1,4) = Dcpg(1,a)
        Hi(2,1) = 1
        Hi(2,2) = 0
        Hi(2,3) = 0
        Hi(2,4) = 0

        THi(1,1) = Hi(1,1)
        THi(2,1) = Hi(1,2)
        THi(3,1) = Hi(1,3)
        THi(4,1) = Hi(1,4)
        THi(1,2) = Hi(2,1)
        THi(2,2) = Hi(2,2)
        THi(3,2) = Hi(2,3)
        THi(4,2) = Hi(2,4)

        RTM = matmul(matmul(Gi,W),TGi)
        dnum = 1
        enum = 1
        do while(dnum<=4)
            enum = 1
            do while(enum<=4)
                Amat(dnum,enum) = -1* F(dnum,enum) * dt
                Amat(4+dnum,enum) = 0 * dt
                Amat(dnum,enum+4) = RTM(dnum,enum) * dt
                Amat(dnum+4,enum+4) = TF(dnum,enum) * dt
                enum = enum+1
            end do
            dnum = dnum+1
        end do

        dnum = 1
        summat = Iden
        val = Iden
        do while (dnum < 300)
            val = (matmul(val, Amat))/dnum
            summat = summat + val
            dnum = dnum + 1
        end do

        dnum = 1
        enum = 1
        do while (dnum<=4)
            enum = 1
            do while(enum<=4)
                RTM(dnum,enum) = summat(dnum,4+enum)
                RTB(dnum,enum) = summat(4+dnum,4+enum)
                enum = enum + 1
            end do
            dnum = dnum + 1
        end do

        dnum = 1
        enum = 1
        do while (dnum<=4)
            enum = 1
            do while(enum<=4)
                TRTB(enum,dnum) = RTB(dnum, enum)
                enum = enum + 1
            end do
            dnum = dnum + 1
        end do

        Q = matmul(TRTB,RTM)

        P = matmul((matmul(F,P)),TF) + Q
        tex = matmul(matmul(Hi,P),THi) + Ri
        lol = inv(tex)
        KF1 = matmul((matmul(P,THi)),lol)

        t1(1,1) = x(1,a)
        t1(2,1) = x(2,a)
        t1(3,1) = x(3,a)
        t1(4,1) = x(4,a)
        temp3(1,1) = zi(1,a)
        temp3(2,1) = zi(2,a)
        t2 = matmul(KF1,(temp3 - matmul(Hi,t1)))

        x(1,a) = t1(1,1) + t2(1,1)
        x(2,a) = t1(2,1) + t2(2,1)
        x(3,a) = t1(3,1) + t2(3,1)
        x(4,a) = t1(4,1) + t2(4,1)
        P = matmul((I-matmul(KF1,Hi)),P)
        a = a + 1
    end do
    a = 1
    do while(a<=Nt)
        Cmmat(b,a) = x(3,a)
        Cdmat(b,a) = x(4,a)
        print*,x(3,a)
        a = a + 1
    end do
    b = b + 1
end do
a = 1
b = 1
sumCd = 0
sumCm = 0
do while(a<=Nt)
    b = 1
    sumCm = 0
    sumCd = 0
    do while(b<=Numi)
        sumCd = sumCd + Cdmat(b,a)
        sumCm = SumCm + Cmmat(b,a)
        b = b+1
    end do
    Cmfinal(1,a) = sumCm/num
    Cdfinal(1,a) = sumCd/num
    a = a+1
end do

!open (unit = 1, file = "Plot_xtsp.dat")
!open (unit = 2, file = "Plot_xtsv.dat")
!open (unit = 3, file = "Plot_IC.dat")
!open (unit = 4, file = "Plot_DC.dat")
open (unit = 5, file = "Plot_Savex.dat")
!open (unit = 6, file = "Plot_Savext2.dat")
!open (unit = 7, file = "Plot_znoise.dat")
open (unit = 8, file = "Plot_CM.dat")
open (unit = 9, file = "Plot_CD.dat")
!open (unit = 10, file = "Plot_Value.dat")

a = 1
do while(a<=Nt)
    !write (1,*) xts(1,a)
    !write (2,*) xts(2,a)
    !write (3,*) Icpg(1,a)
    !write (4,*) Dcpg(1,a)
    write (5,*) x(1,a)
    !write (6,*) xt(2,a)
    !write (7,*) zi(1,a)
    write (8,*) Cmfinal(1,a)
    write (9,*) Cdfinal(1,a)
    !write (10,*) value(1,a)
    a = a+1
end do

END PROGRAM
