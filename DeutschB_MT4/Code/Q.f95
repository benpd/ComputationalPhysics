module mod

	use numtype
	use integr
	implicit none
	
	real(dp), parameter :: w    = 1._dp
	real(dp), parameter :: hbar = 1._dp
	real(dp), parameter :: mass = 1._dp
 	integer :: i, k, j, m, info
	integer, parameter :: ndim = 10, lwork= 16*ndim 
	real(dp) :: a, b, c, d, dia, work(lwork) 
	real(dp) :: res, resdos
	integer :: nint, ifail
	real(dp), dimension(maxint) :: xint, wint

	real(dp), dimension (0:ndim, 0:ndim) :: E, BigU, BigV, T, H
	real(dp), dimension (1:ndim+1, 1:ndim+1) :: HH
	real(dp), dimension (1:ndim+1) :: O

end module mod

Program BUharmo
	
	use mod
	real(dp), external :: psiguy2
	real(dp), external :: U, V

!MAKING EXPLICIT E MATRIX
		do  m = 0, ndim
			do  j = 0, ndim 
				dia = 0.5_dp
				if ( m==j ) then                                                                                 
                	E(m,j) = (dia + j)*w
               
            	else 
					E(m,j) = 0._dp
            
				end if 

			end do 
		end do
	! do i =0, ndim
	! 	print '(''(''11(f8.3, 3x)'')'')', E(i, 0:ndim)
	! end do 

!INTEGRATION
	nint = 100
	a = -10._dp
    b = 10._dp
	c = 0._dp
	d = 0._dp

    call d01bcf(0,a,b,c,d,nint,wint,xint,ifail)

	do i = 0, ndim
		do j= 0, ndim 
			res = 0._dp
			resdos = 0._dp
			do k = 0, nint
        		res = res + wint(k)*(psiguy2(xint(k),i)) * psiguy2(xint(k), j) * U(xint(k))  
				resdos = resdos + wint(k)*(psiguy2(xint(k),i)) * psiguy2(xint(k), j) * V(xint(k))  
			end do
			bigU(i, j) = res
			bigV(i, j) = resdos
		end do
	end do 	

!CREATION OF BIG ELEMENTS OF EQUATIONS
print *, '_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'	
 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', bigU(i, 0:ndim)
	end do 

print *, '_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'	
 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', bigV(i, 0:ndim)
	end do 

!CREATION OF MATRIX EQUATIONS
T(0:ndim, 0:ndim)= E(0:ndim, 0:ndim)- bigU(0:ndim, 0:ndim)
H(0:ndim, 0:ndim)= T(0:ndim, 0:ndim)+ bigV(0:ndim, 0:ndim)
HH(1:ndim+1, 1:ndim+1) = H(0:ndim, 0:ndim)

 print *, '_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'	

 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', H(i, 0:ndim)
	end do 

call dsyev('v', 'u', ndim+1, HH, ndim+1, O, work, lwork, info)

print *, '_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'
	do i = 1, ndim
		if(O(i)<0) then	
			print *, i, O(i)
		end if	
	end do

print *, '_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'	
 	do i =1, ndim+1
		print '(''(''11(f8.3, 3x)'')'')', HH(i, 1:ndim+1)
	end do

End program BUharmo

!EXTERNAL FUNCTIONS 
Function psiguy2(x,n)
	use mod 
 	integer :: n,l
	real(dp) :: x, psiguy2
	real(dp), dimension(0:n) :: psi
	
		psi(0) = (mass*w/(hbar*pi))**(0.25_dp) * exp(-mass*w*x**2/(2*hbar))
		
	do l = 0, n-1	
 		psi(l+1) = (2*sqrt(mass*w/hbar)*x)/(2._dp*(l+1))**(0.5_dp)*psi(l)& 
 				 		   -(float(l)/(l+1))**(0.5_dp)*psi(l-1) 	 
	end do
			
		psiguy2 = psi(n)

End function 

Function U(x) 
	use mod 
	real(dp) :: U, x

		U= 0.5_dp * mass * (w**2) * x**2 

End function 

Function V(x) 
	use mod
	real(dp) :: V, x

		V= -10._dp/abs(iic+x)

End function 