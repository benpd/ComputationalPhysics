Program BB
	
	use numtype 
	implicit none
	
	real(dp), parameter ::  c = 3.00E8 		  	! speed of light
	real(dp), parameter ::	h = 6.62E-34       	! Planck's constant  
	real(dp), parameter ::  kb = 1.38E-23    	! Boltzmann constant

	real(dp) :: T1 = 2.7						! temps	in Kelvin
	real(dp) :: T2 = 300						 
	real(dp) :: T3 = 6500						 
	
	real(dp) :: v   							! frequency absorbed in EE
	real(dp) :: constant = 2* (kb**3)/(c*h)**2	! constant package 
	real(dp) :: EE 								! new replace for E 
	real(dp) :: Eq 							 	! result from equation 1   
	
	integer:: i									! counting integer 
	
	
!_________TEMP 1___________		
	do i = 0, 10000
	EE  = float(i) * 0.001 
	Eq = constant * T1**3 * EE**3 / (exp(EE) - 1) 
	
	v   = EE * kb * T1/h 
	
	print *, EE, Eq
	write (11,*) EE, Eq 
	write (21,*) v, Eq		
	
	Eq = Eq * 3E10								! scale factor to 33 
	Write (31,*) EE, Eq
	
	end do 
	
	
!_________TEMP 2___________		
	do i = 0, 10000
	EE  = float(i) * 0.001 
	Eq = constant * T2**3 * EE**3 / (exp(EE) - 1) 
	
	v   = EE * kb * T1/h 
	
	print *, EE, Eq 
	write (12,*) EE, Eq 
	write (22,*) v, Eq
	
	Eq = Eq * 1.5E4								!scale factor to 33
	Write (32,*) EE, Eq
	
	end do 
	
	
!_________TEMP 3___________	
	do i = 0, 10000
	EE  = float(i) * 0.001 
	Eq = constant * T3**3 * EE**3 / (exp(EE) - 1) 
	
	v   = EE * kb * T1/h 
	
	print *, EE, Eq 
	write (13,*) EE, Eq							! NO SCALE FACTOR 
	write (23,*) v, Eq 							
	
	
	end do 


end program BB