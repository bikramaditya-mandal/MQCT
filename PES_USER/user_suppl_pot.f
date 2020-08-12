!----------!
! OPTION 1 !
!----------! 
!	  THIS IS THE DEFAULT OPTION IN THE CODE, NO KEYWORDS ARE REQUIRED 
      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT NONE
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
      STOP "ERROR: USER_DEFINED_PES IS NOT SUPPLIED"
      END SUBROUTINE USER_DEFINED_PES
!----------!
! OPTION 2 !
!----------! 
!	  USE KEYWORD "EXPANSION=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_TERMS(T,I,R)
!     THIS SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION AT A GIVEN DISTANCE
!     INPUT:  R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent 	  
      IMPLICIT NONE	  
      REAL*8 T,R
      INTEGER I
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
      STOP "ERROR: USER_DEFINED_TERMS IS NOT SUPPLIED"
	  END SUBROUTINE USER_DEFINED_TERMS
!----------!
! OPTION 3 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_FILE=YES" TO INITIATE THIS OPTION
! 	  SIMILAR TO OPTION 2, BUT NO SUBROUTINE IS REQUIRED
!     USER SHOULD PROVIDE THE FILE EXPAN_PES_TERMS.DAT 
!     IN THE MAIN PROGRAM DIRECTORY CONTAINING THE COEFFICEINS 
!     OF POTENTIAL EXPANSION PRECOMPUTED EXTERNALLY.
! 	  SEE EXAMPLE FILES SUPPLIED WITH THE CODE.
!----------!
! OPTION 4 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_ONFLY=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_COEFFS(T,DTDR,I,R) 
!     THIS SUBROUTINE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION 
!     AND THEIR DERIVATIVES AT A GIVEN DISTANCE R
      IMPLICIT NONE
!     INPUT : R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent, DTDR - its radial derivative 	  
      REAL*8 T,R,DTDR 
      INTEGER I
!     USER MUST INCERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:	
      STOP "ERROR: USER_DEFINED_COEFFS IS NOT SUPPLIED"
      END SUBROUTINE USER_DEFINED_COEFFS 
!----------! 
		subroutine MFTOBF_CONV(a1,b1,g1,a2,b2,g2,theta,phi,thetap,phip)
		implicit none
		real*8 o(3),h1(3),h2(3)
		real*8 safe_o(3),safe_h1(3),safe_h2(3)
		real*8 h3(3),h4(3)
		real*8 com_h2(3),r_com_h2,r_com_h,temp_com_h2(3)
		real*8 tempo(3),temph1(3),temph2(3)
		real*8 temph3(3),temph4(3)
		real*8 r_oh1,r_oh2,r_h1h2
		real*8 rr_oh1,rr_oh2,rr_h1h2
		real*8 temp1(3,3),temp2(3),temp3(3)
		real*8 theta,phi,temp00,temp01,temp02,temp03,temp04
		real*8 thetap,phip,temp10,temp11,temp12,temp13,temp14
		real*8,parameter :: r_w_h=6.00d0
		real*8,parameter :: r_h2=0.740d0
		real*8,parameter :: pi=4.0d0*datan(1.0d0)
		real*8 a1,b1,g1,rot_mat1(3,3),inv_rot_mat1(3,3)
		real*8 temp_rot_mat1(3,3)
		real*8 a2,b2,g2,rot_mat2(3,3)
		
		! coordinates of H2O
		o(1)=0.0d0
		o(2)=0.0d0
		o(3)=0.06530d0
		h1(1)=0.774880d0
		h1(2)=0.0d0
		h1(3)=-0.53470d0
		h2(1)=-0.774880d0
		h2(2)=0.0d0
		h2(3)=-0.53470d0
		
		safe_o=o
		safe_h1=h1
		safe_h2=h2
		
		! coordinates of hydrogens in H2
		h3(1)=0.0d0
		h3(2)=0.0d0
		h3(3)=r_h2*0.50d0
		
		h4(1)=0.0d0
		h4(2)=0.0d0
		h4(3)=-r_h2*0.50d0
		
		! find the COM of H2 molecule
		com_h2(1)=(h3(1)+h4(1))*0.50d0
		com_h2(2)=(h3(2)+h4(2))*0.50d0
		com_h2(3)=(h3(3)+h4(3))*0.50d0

		!find inter-atomic distance of H2O
		rr_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		rr_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		rr_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		!Euler rotation matrix for H2O
		call rot_mat(a1,b1,g1,rot_mat1)
		
		!Euler rotation matrix for H2
		call rot_mat(a2,b2,g2,rot_mat2)
		
		!SF rotation matrix for H2
		call sf_rot_mat(b2,a2,rot_mat2)
		
		tempo=o
		temph1=h1
		temph2=h2
		o=matmul(rot_mat1,tempo)
		h1=matmul(rot_mat1,temph1)
		h2=matmul(rot_mat1,temph2)
		
		!find inter-atomic distance of H2O after rotation
		r_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		r_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		r_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		if(abs(rr_oh1-r_oh1).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & O_H1",rr_oh1,r_oh1
		stop
		endif
		if(abs(rr_oh2-r_oh2).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & O_H2",rr_oh2,r_oh2
		stop
		endif
		if(abs(rr_h1h2-r_h1h2).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & H1_H2",rr_h1h2,r_h1h2
		stop
		endif

		call inv_rot_mat(a1,b1,g1,inv_rot_mat1)
	!	call inverse(rot_mat1,inv_rot_mat1,3)
	!	call rot_mat(-a1,-b1,-g1,inv_rot_mat1)
		tempo=o
		temph1=h1
		temph2=h2
		o=matmul(inv_rot_mat1,tempo)
		h1=matmul(inv_rot_mat1,temph1)
		h2=matmul(inv_rot_mat1,temph2)
		
		!find inter-atomic distance of H2O after inverse-rotation
		r_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		r_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		r_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		
		if(abs(rr_oh1-r_oh1).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & O_H1",rr_oh1,r_oh1
		stop
		endif
		if(abs(rr_oh2-r_oh2).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & O_H2",rr_oh2,r_oh2
		stop
		endif
		if(abs(rr_h1h2-r_h1h2).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & H1_H2",rr_h1h2,r_h1h2
		stop
		endif
		
		!check if H2O is shifted back to its reference position
		if(abs(safe_o(1)-o(1)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & X-coordinate",safe_o(1),o(1)
		stop
		endif
		if(abs(safe_o(2)-o(2)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & Y-coordinate",safe_o(2),o(2)
		stop
		endif
		if(abs(safe_o(3)-o(3)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & Z-coordinate",safe_o(3),o(3)
		stop
		endif
		if(abs(safe_h1(1)-h1(1)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & X-coordinate",safe_h1(1),h1(1)
		stop
		endif
		if(abs(safe_h1(2)-h1(2)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & Y-coordinate",safe_h1(2),h1(2)
		stop
		endif
		if(abs(safe_h1(3)-h1(3)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & Z-coordinate",safe_h1(3),h1(3)
		stop
		endif
		if(abs(safe_h2(1)-h2(1)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & X-coordinate",safe_h2(1),h2(1)
		stop
		endif
		if(abs(safe_h2(2)-h2(2)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & Y-coordinate",safe_h2(2),h2(2)
		stop
		endif
		if(abs(safe_h2(3)-h2(3)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & Z-coordinate",safe_h2(3),h2(3)
		stop
		endif
		
		! rotation on H2 molecule
		temph3=h3
		temph4=h4
		
		h3=matmul(rot_mat2,temph3)
		h4=matmul(rot_mat2,temph4)
		
		! shifting of hydrogens of H2
		h3(3)=h3(3)+r_w_h
		h4(3)=h4(3)+r_w_h	
		com_h2(3)=com_h2(3)+r_w_h
		
		temph3=h3
		temph4=h4
		temp_com_h2=com_h2
		h3=matmul(inv_rot_mat1,temph3)
		h4=matmul(inv_rot_mat1,temph4)	
		com_h2=matmul(inv_rot_mat1,temp_com_h2)	
		
		r_com_h2=dsqrt(com_h2(1)**2+com_h2(2)**2+com_h2(3)**2)
		temp01=com_h2(3)/r_com_h2
		theta=dacos(temp01)
		
		phi=datan2(com_h2(2),com_h2(1))
		if(phi.lt.0.00d0) phi=phi+2.0d0*pi
		
		if(theta.ne.theta) then
		print*,'Something is wrong in the transformation 
     & of angle theta',theta
		stop
		endif
		if(phi.ne.phi) then
		print*,'Something is wrong in the transformation 
     & of angle phi',phi
		stop
		endif
		
		! shifting COM of H2 and H2 molecule to origin
		h3=h3-com_h2
		h4=h4-com_h2
		
		r_com_h=dsqrt(h3(1)**2+h3(2)**2+h3(3)**2)
		temp11=h3(3)/r_com_h
		thetap=dacos(temp11)

		phip=datan2(h3(2),h3(1))
		if(phip.lt.0.00d0) phip=phip+2.0d0*pi
		
		if(thetap.ne.thetap) then
		print*,'Something is wrong in the transformation 
     & of angle theta_prime',thetap
		stop
		endif
		if(phip.ne.phip) then
		print*,'Something is wrong in the transformation 
     & of angle phi_prime',phip
		stop
		endif
		
		contains
		
		subroutine rot_mat(a,b,g,mat)
		implicit none
		real*8 a,b,g,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(a)*dcos(b)*dcos(g)-dsin(a)*dsin(g)
		mat(1,2)=-dcos(a)*dcos(b)*dsin(g)-dsin(a)*dcos(g)
		mat(1,3)=dcos(a)*dsin(b)
		
		mat(2,1)=dsin(a)*dcos(b)*dcos(g)+dcos(a)*dsin(g)
		mat(2,2)=-dsin(a)*dcos(b)*dsin(g)+dcos(a)*dcos(g)
		mat(2,3)=dsin(a)*dsin(b)
		
		mat(3,1)=-dsin(b)*dcos(g)
		mat(3,2)=dsin(b)*dsin(g)
		mat(3,3)=dcos(b)
		endsubroutine
		
		subroutine inv_rot_mat(a,b,g,mat)
		implicit none
		real*8 a,b,g,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(a)*dcos(b)*dcos(g)-dsin(a)*dsin(g)
		mat(2,1)=-dcos(a)*dcos(b)*dsin(g)-dsin(a)*dcos(g)
		mat(3,1)=dcos(a)*dsin(b)
		
		mat(1,2)=dsin(a)*dcos(b)*dcos(g)+dcos(a)*dsin(g)
		mat(2,2)=-dsin(a)*dcos(b)*dsin(g)+dcos(a)*dcos(g)
		mat(3,2)=dsin(a)*dsin(b)
		
		mat(1,3)=-dsin(b)*dcos(g)
		mat(2,3)=dsin(b)*dsin(g)
		mat(3,3)=dcos(b)
		endsubroutine

		subroutine sf_rot_mat(sf_t,sf_p,mat)
		implicit none
		real*8 sf_t,sf_p,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(sf_p)*dcos(sf_t)
		mat(1,2)=-dsin(sf_p)
		mat(1,3)=dcos(sf_p)*dsin(sf_t)
		
		mat(2,1)=dsin(sf_p)*dcos(sf_t)
		mat(2,2)=dcos(sf_p)
		mat(2,3)=dsin(sf_p)*dsin(sf_t)
		
		mat(3,1)=-dsin(sf_t)
		mat(3,2)=0.00d0
		mat(3,3)=dcos(sf_t)
		endsubroutine
		
		subroutine mat_vec(row,column,mat,vec,res)
		implicit none
		integer i,j
		integer row,column
		real*8 mat(row,column),vec(column),res(row)
		
		res=0.00d0
		do i=1,row
		res(i)=0.00d0
		do j=1,column
		res(i)=res(i)+mat(i,j)*vec(j)
		enddo
		enddo
		
		end subroutine
		
		subroutine mat_print(m,n,a)
		implicit none
		integer m,n,i,j
		real*8 a(m,n)
		
		do i=1,m
		do j=1,n
		if(abs(a(i,j)).le.1d-9 .and. a(i,j).ge.0.00d0) a(i,j)=0.00d0
		enddo
		enddo
		
		do i=1,m
		do j=1,n
		write(*,'(f12.5,2x)',advance='no')a(i,j)
		enddo
		write(*,*)' '
		enddo
		write(*,*)' '
		endsubroutine
		
		subroutine inverse(a,c,n)
		!============================================================
		! Inverse matrix
		! Method: Based on Doolittle LU factorization for Ax=b
		! Alex G. December 2009
		!-----------------------------------------------------------
		! input ...
		! a(n,n) - array of coefficients for matrix A
		! n      - dimension
		! output ...
		! c(n,n) - inverse matrix of A
		! comments ...
		! the original matrix a(n,n) will be destroyed 
		! during the calculation
		!===========================================================
		implicit none 
		integer n
		double precision a(n,n), c(n,n)
		double precision L(n,n), U(n,n), b(n), d(n), x(n)
		double precision coeff
		integer i, j, k

		! step 0: initialization for matrices L and U and b
		! Fortran 90/95 aloows such operations on matrices
		L=0.0
		U=0.0
		b=0.0

		! step 1: forward elimination
		do k=1, n-1
		   do i=k+1,n
			  coeff=a(i,k)/a(k,k)
			  L(i,k) = coeff
			  do j=k+1,n
				 a(i,j) = a(i,j)-coeff*a(k,j)
			  end do
		   end do
		end do

		! Step 2: prepare L and U matrices 
		! L matrix is a matrix of the elimination coefficient
		! + the diagonal elements are 1.0
		do i=1,n
		  L(i,i) = 1.0
		end do
		! U matrix is the upper triangular part of A
		do j=1,n
		  do i=1,j
			U(i,j) = a(i,j)
		  end do
		end do

		! Step 3: compute columns of the inverse matrix C
		do k=1,n
		  b(k)=1.0
		  d(1) = b(1)
		! Step 3a: Solve Ld=b using the forward substitution
		  do i=2,n
			d(i)=b(i)
			do j=1,i-1
			  d(i) = d(i) - L(i,j)*d(j)
			end do
		  end do
		! Step 3b: Solve Ux=d using the back substitution
		  x(n)=d(n)/U(n,n)
		  do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
			  x(i)=x(i)-U(i,j)*x(j)
			end do
			x(i) = x(i)/u(i,i)
		  end do
		! Step 3c: fill the solutions x(n) into column k of C
		  do i=1,n
			c(i,k) = x(i)
		  end do
		  b(k)=0.0
		end do
		end subroutine inverse 
		
		end subroutine
		