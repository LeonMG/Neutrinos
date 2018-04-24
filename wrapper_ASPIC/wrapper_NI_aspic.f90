! =======================================================================================
!
!		Wrapper that calculates the primordial power spectrum in the
!		     Natural Inflation model (NI) from the ASPIC library
!
!				     Version 1 - March 2018
!			    Gabriela Aguilar (gaguilar@astro.unam.mx)
!				  Instituto de Astronomia, UNAM
!
!
!    The input parameters are:
!	A_s    = amplitude of scalar fluctuations
!	lnRrad = reheating parameter 
!	log10f = model parameter for natural inflation ( = log10(f/Mpl) )
!
!
!    The columns of the output file are:
!	1st column = k in Mpc^-1 (log-spaced from 1e-6 to 1)
!	2nd column = primordial scalar power spectrum
!	3rd column = primordial tensor power spectrum
!
!
! =======================================================================================


 program main


     include "aspic.h"

     use infprec, only : kp
     use cosmopar, only : lnRhoNuc, powerAmpScalar
     use nisr, only : ni_epsilon_one, ni_epsilon_two, ni_epsilon_three
     use nireheat, only : ni_x_rrad
 
     implicit none


     ! Some constants for the code
     Integer,  parameter :: npk = 20				! Number of points for the log-spaced of k
!     Real(kp), parameter :: PI = 4._kp*DATAN(1._kp)		! pi
!     Real(kp), parameter :: Ggrav = 				! Gravitational constant
!     Real(kp), parameter :: Mpl2 = ( 8._kp * PI * Ggrav )**2	! Squared reduced Planck mass
     Real(kp), parameter :: Cinai = - 0.7296_kp			! Cinai = gamma_E + ln(2) - 2 =(aprox.) -0.7296, gamma_E = euler constant
     Real(kp), parameter :: kpivot = 0.05_kp			! Mpc^-1
     Real(kp), parameter :: kmin = 1.e-6_kp
     Real(kp), parameter :: kmax = 1._kp


     Integer :: i
     Real(kp) :: A_s, lnRrad, log10f				! A_s = amplitude of scalar fluctuations, InRad = reheating parameter, logf =  model parameter for Natural Inflation
     Real(kp) :: xstar, f					! xstar = phi/Mpl, f = f/Mpl
     Real(kp) :: eps1, eps2, eps3				! epsilons in eqs. 2.20 to 2.25 of 1303.3787
     Real(kp) :: Pstar						! Pstar = amplitude of the scalar power spectrum evaluated at 'kstar', the pivot wave-number (default 0.05 Mpc^-1)
     Real(kp) :: a0s, a1s, a2s, a0T, a1T, a2T 
     Real(kp) :: C2inai, ginai, finai, PI2
     Real(kp) :: kwn, lk, lkmin, lkmax, dlk
     Real(kp) :: P_s, P_T, P0_s, P0_T
     Character(len=20) :: filename


     ! 1. Obtain the input parameters
     Namelist /Input/ A_s, lnRrad, log10f


     ! 2. Read the input parameters from the input file
     open (3, file='iii.par', status = 'old' )
     read (3, nml = Input)
     close(3)


     ! 3. Name and open file for the output file
     filename = 'output_PS.dat'
     open(10,file=filename,form='formatted',status='replace')


     ! 4. Give some parameters for the functions
     Pstar  = powerAmpScalar
     f      = 10._kp**(log10f)				! log10(f/Mpl) = log10f -> f/Mpl = 10^(log10f)
     ginai  = 7._kp
     finai  = 5._kp
     C2inai = Cinai * Cinai
     PI2    = PI * PI
     lkmin  = log10(kmin)
     lkmax  = log10(kmax)
     dlk    = (lkmax - lkmin) / real(npk-1,kp)


     ! 5. Calculate  and epsilons in eqs. 2.20 t0 2.25 of 1303.3787
     xstar = ni_x_rrad(f,lnRrad,Pstar)
     eps1  = ni_epsilon_one(xstar,f)			! = eps1/Mpl^2
     eps2  = ni_epsilon_two(xstar,f)
     eps3  = ni_epsilon_three(xstar,f)

     
     ! 6. Calculate the expansion parameters in 2.18 of 1303.3787, for the scalar power spectrum
     a0s = 1._kp - 2._kp * ( Cinai + 1._kp ) * eps1 - Cinai * eps2 &		
         + ( 2._kp*C2inai + 2._kp*Cinai + 0.5_kp*PI2 - finai ) * eps1**2 &
	 + ( C2inai - Cinai + 0.583333_kp*PI2 - ginai ) * eps1 * eps2 &
	 + ( 0.5_kp*C2inai + 0.125_kp*PI2 - 1._kp ) * eps2**2 &
	 + ( - 0.5_kp*C2inai + PI2/24._kp ) * eps2 * eps3				! eq. 2.20 of 1303.3787

     a1s = - 2._kp * eps1 - eps2 + 2._kp * ( 2._kp*Cinai + 1._kp ) * eps1**2 &
	 + ( 2._kp*Cinai - 1._kp ) * eps1*eps2 + Cinai*eps2**2 - Cinai*eps2*eps3	! eq. 2.21 of 1303.3787
 
     a2s = 4._kp*eps1**2 + 2._kp*eps1*eps2 + eps2**2 - eps2*eps3			! eq. 2.22 of 1303.3787


     ! 7. Calculate the expansion parameters in 2.18 of 1303.3787, for the tensor power spectrum
     a0T = 1._kp - 2._kp * ( Cinai + 1._kp ) * eps1 &
	 + ( 2._kp*C2inai + 2._kp*Cinai + 0.5_kp*PI2 - finai ) * eps1**2 &
	 + ( - C2inai - 2._kp*Cinai + PI2 / 12._kp - 2._kp ) * eps1 * eps2		! eq. 2.23 of 1303.3787

     a1T = - 2._kp * eps1 + 2._kp * ( 2._kp*Cinai + 1._kp ) * eps1**2 &
	 - 2._kp * ( Cinai + 1._kp ) * eps1 * eps2					! eq. 2.24 of 1303.3787

     a2T = 4._kp * eps1**2 - 2._kp * eps1 * eps2					! eq. 2.25 of 1303.3787


     ! 8. Calculate the P0 quantity in 2.18 of 1303.3787, for scalar and tensor power spectrum, respectively
     P0_s = A_s / a0s					! scalar
     P0_T = 16._kp * eps1 * P0_s			! Tensor


     ! 9. Compute the scalar and tensor power spectra with eq. 2.18 of 1303.3787 
     do i=1, npk

        ! 9.1 log-spaced wave number
	lk  = lkmin + real(i-1,kp) * dlk
	kwn = 10._kp**lk

	! 9.2 Scalar power spectrum
	call power_spectrum(P_s,kwn,P0_s,a0s,a1s,a2s,kpivot)

	! 9.3 Tensor power spectrum
	call power_spectrum(P_T,kwn,P0_T,a0T,a1T,a2T,kpivot)

	! 9.4 Write the colums of the output file
	write(10,*) kwn, P_s, P_T

     end do


     ! 10. Close the output file
     close(10)


     ! End of the program

 end program main





! -------------------------------------------------------------------------------------------------
! -------------------------   Auxiliary subroutines for the main program   ------------------------
! -------------------------------------------------------------------------------------------------

 subroutine power_spectrum(PSk,kwni,P0,a0,a1,a2,kpiv)

     ! This routine calculates the Power spectrum using eq. 2.18 of 1303.3787

     include "aspic.h"

     use infprec, only : kp


     Real(kp), intent(in) :: kwni, P0, a0, a1, a2, kpiv
     Real(kp), intent(inout) :: Psk

     Real(kp) :: aux


     ! Auxiliary quantity
     aux = log(kwni/kpiv)


     ! Power spectrum from eq. 2.18 of 1303.3787
     Psk = P0 * ( a0 + a1 * aux + 0.5_kp * a2 * aux**2 )


 end subroutine power_spectrum






