
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preconditionings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	subroutine getPrecondPara(VM2,AM2,RR,TT,RMIUT,ALAGMR,ubetac2)
	
	Use Global
	IMPLICIT NONE

	REAL:: VM2,AM2,RR,TT,RMIUT,ALAGMR
	REAL:: ubetac2
	
	REAL:: K,Kc,RmiuL,RmiuAll

		K=3.0 !2.0_8
		Kc=1.0 !4.0_8/3.0_8

		ubetac2 = min(max(VM2/AM2, K*XM*XM),1.0)
!			ubetac2 = min(vm2/am2, 1.0_8) ! (Ur/c)**2 c=am
!			ubetac2 = max(ubetac2, 1.0d-10) ! !1.0d-5**2   !or max(ubetac2, Xm2)???

!		viscosity limit
		 !TTL= am2/ga/PPF !Tm
         RmiuL=(1.0_8+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)
		 RmiuAll = RmiuL+ RmiuT !VL(1,8,I,J,K)
		 ubetac2=max(ubetac2, Kc*(RmiuALL/RR/Ref/ALAGMR)**2.0/AM2)

	return
	end subroutine




	

	subroutine DWtoDQtoDW (AinRo,AinRx,AinRy,AinRz,AinRe, &
     &			pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT,  &
     &			AinP,AinX,AinY,AinZ,AinT)

	Use Global
	IMPLICIT NONE

	!include "head\000_head_common.fd"
	real:: AinRo,AinRx,AinRy,AinRz,AinRe
	real:: pRo,pVx,pVy,pVz,pE, pPP,pTT
	real:: AAc,ubetac2
	real:: AinP,AinX,AinY,AinZ,AinT
	
	REAL:: roinv,rCp,a11,a15,b11,b15

		!K transform
			AinRe = AinRe - (pVx*AinRx+pVy*AinRy+pVz*AinRz)     &
     &			-(pE/1.4-1.2*(pVx*pVx+pVy*pVy+pVz*pVz))*AinRo
			AinRx = AinRx - pVx*AinRo
			AinRy = AinRy - pVy*AinRo
			AinRz = AinRz - pVz*AinRo
			AinRo = AinRo

		!beta transform
			roinv= 1.0/pRo

			rCp = 1.4/0.4*PPF * pRo
			!a_det = rCp/(AAc*AAc*ubetac2)   !=rCp/a11

			a11=AAc*AAc*ubetac2
			a15=-0.4*ubetac2/rCp  !-0.4/AAc/AAc/a_det
			b11=1.0*a11/rCp !1.0/a_det
			b15=(1.0+0.4*ubetac2)/rCp

			AinP=a11*AinRo+a15*AinRe
			AinX=AinRx*roinv
			AinY=AinRy*roinv
			AinZ=AinRz*roinv
			AinT=b11*AinRo+b15*AinRe

			AinRo = AinP/PPF/pTT-pRo*AinT/pTT
			AinRx = AinRo*pVx+pRo*AinX
			AinRy = AinRo*pVy+pRo*AinY
			AinRz = AinRo*pVz+pRo*AinZ
			AinRe = AinP/0.4+(pVx*AinRx+pVy*AinRy+pVz*AinRz)

	return
	end subroutine





	subroutine PrecondDW (AinRo,AinRx,AinRy,AinRz,AinRe,  &
     &			pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT)   !,!     >			AinP,AinX,AinY,AinZ,AinT)

	Use Global
	IMPLICIT NONE

	real AinRo,AinRx,AinRy,AinRz,AinRe
	real pRo,pVx,pVy,pVz,pE
	real AAc,ubetac2
	real pPP,pTT
	real AinP,AinX,AinY,AinZ,AinT
	real tDP, pU2, pUAU, tCDP
	REAL:: rCp,a11,a15,b11,b15
	
		!K transform
		pU2=pVx*pVx+pVy*pVy+pVz*pVz
		pUAU=pVx*AinRx+pVy*AinRy+pVz*AinRz

			AinRe = AinRe - (pUAU)-(pE/1.4-1.2*(pU2))*AinRo

		!beta transform

			rCp = 1.4/0.4*PPF * pRo
			!a_det = rCp/(AAc*AAc*ubetac2)   !=rCp/a11

			a11=AAc*AAc*ubetac2
			a15=-0.4*ubetac2/rCp  !-0.4/AAc/AAc/a_det
			b11=1.0*a11/rCp !1.0/a_det
			b15=(1.0+0.4*ubetac2)/rCp

			AinP=a11*AinRo+a15*AinRe
			AinT=b11*AinRo+b15*AinRe

			AinRo = AinP/PPF/pTT-pRo*AinT/pTT
			AinRe = AinP/0.4+(pUAU)


	return
	end subroutine



!	subroutine PrecondDW3 (AinRo,AinRx,AinRy,AinRz,AinRe,
!     >			pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!!     >			AinP,AinX,AinY,AinZ,AinT)
!!	include "head\000_head_common.fd"
!	real*8 AinRo,AinRx,AinRy,AinRz,AinRe
!	real*8 pRo,pVx,pVy,pVz,pE
!	real*8 AAc,ubetac2
!	real*8 pPP,pTT
!!	real*8 AinP,AinX,AinY,AinZ,AinT
!	real*8 tDP, pU2, tCDP

!	pU2=pVx*pVx+pVy*pVy+pVz*pVz
!	tDP=0.4_8*(0.5_8*pU2* AinRo
!     >			- pVx*AinRx - pVy*AinRy - pVz*AinRz
!     >			+ AinRe)

!	tCDP=(1.0_8-ubetac2)/AAc/AAc*tDP

!	AinRo=AinRo- tCDP
!	AinRx=AinRx- tCDP * pVx
!	AinRy=AinRy- tCDP * pVy
!	AinRz=AinRz- tCDP * pVz
!	AinRe=AinRe- tCDP * (pE+0.5_8*pU2)


!	return
!	end subroutine
