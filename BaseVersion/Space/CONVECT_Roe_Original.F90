SUBROUTINE CONVECT_Roe_Original
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: KKtm
    
    REAL:: CFL_max, CFL_min
    REAL:: dT_max, dT_min
    REAL:: Den,Vxi,Vyi,Vzi
    REAL:: AA,UU,VV,WW,UU1,VV1,WW1
    REAL:: ST1,ST11,ST2,ST22,ST3,ST33
    REAL:: Dtm_local, CFL_local

    REAL:: Gradient
    REAL:: SD1,SD2,SD3,SD4
    REAL:: Kx,Ky,Kz,Kt,Kapa
    REAL:: dn02,dn01,dn00,dn10
    REAL:: Q02,Q01,Q00,Q10
    REAL:: Fl(7),Fr(7),Diff(7)
    REAL:: RRl,RRr,PPl,PPr,UUl,UUr,VVl,VVr,WWl,WWr,rKl,rKr
    REAL:: U00l,U00r,Ub_l,Ub_r
    REAL:: Dr,Du,Dv,Dw,Dp,DU_b
    REAL:: Rroe,Uroe,Vroe,Wroe,Proe,AAroe,Ubroe, Hroe
    REAL:: rKro
    REAL:: ALF1,ALF2,ALF3,ALF4,ALF5,ALF6,ALF7,ALF8

    REAL:: ULo,URo,VLo,VRo,WLo,WRo,PLo,PRo,TLo,TRo
    REAL:: rcpcv,DENLo,DENRo
    REAL:: ga,ga1
    
    REAL:: sav1,sav2,sav3,sav
    REAL:: RmiuT,AlagmR
    REAL:: DroR,DrxR,DryR,DrzR,DreR
    REAL:: TTL
    REAL:: pGam
    real::rLo,rRo
    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1



!!!!!!
!	  Roe scheme with R-S entropy correction
!!!!!!

!	III=IF_PRECONDITION
!	IF_PRECONDITION=0

!     i direction

      do k=1,NK1
      do j=1,NJ1
      do i=1,NI
!        VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm

         ULo=VL(1,1,i,j,k)
	   URo=VR(1,1,i,j,k)
         !!!!!!!!!!!!!!!!!
         VLo=VL(1,2,i,j,k)
	   VRo=VR(1,2,i,j,k)
         !!!!!!!!!!!!!!!!!
         WLo=VL(1,3,i,j,k)
	   WRo=VR(1,3,i,j,k)

         PLo=VL(1,4,i,j,k)
	   PRo=VR(1,4,i,j,k)

	   TLo=VL(1,5,i,j,k)
	   TRo=VR(1,5,i,j,k)
         !!!!!!!!!!!!!!!!!
        rLo=rccl(1,i,j,k)*omega(1)   !by ydd
        rRo=rccr(1,i,j,k)*omega(1)


          sav1=SD(1,1,I,J,K)  !aix(i,j,k)
          sav2=SD(1,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(1,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(1,I,J,K)   !aim(i,j,k)+small
        RmiuT = VL(1,8,I,J,K)
        AlagmR = ALAGM(1,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I-1,J,K))


        CALL CONVECT_Roe_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR,rLo,rRo)  !output


 
	    D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I-1,J,K)=D(1,I-1,J,K)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I-1,J,K)=D(2,I-1,J,K)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I-1,J,K)=D(3,I-1,J,K)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I-1,J,K)=D(4,I-1,J,K)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I-1,J,K)=D(5,I-1,J,K)+dreR * sav

	end do
	end do
	end do


!     j direction

      do k=1,NK1
      do j=1,NJ
      do i=1,NI1

         ULo=VL(2,1,i,j,k)
	   URo=VR(2,1,i,j,k)
         !!!!!!!!!!!!!!!!!
         VLo=VL(2,2,i,j,k)
	   VRo=VR(2,2,i,j,k)
         !!!!!!!!!!!!!!!!!
         WLo=VL(2,3,i,j,k)
	   WRo=VR(2,3,i,j,k)

         PLo=VL(2,4,i,j,k)
	   PRo=VR(2,4,i,j,k)

	   TLo=VL(2,5,i,j,k)
	   TRo=VR(2,5,i,j,k)
         !!!!!!!!!!!!!!!!!
        rLo=rccl(2,i,j,k)*omega(1)
        rRo=rccr(2,i,j,k)*omega(1)


          sav1=SD(2,1,I,J,K)  !aix(i,j,k)
          sav2=SD(2,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(2,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(2,I,J,K)   !aim(i,j,k)+small
		 RmiuT = VL(2,8,I,J,K)
        AlagmR = ALAGM(2,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I,J-1,K))


        CALL CONVECT_Roe_Common( pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR,rLo,rRo)  !output

	    D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I,J-1,K)=D(1,I,J-1,K)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I,J-1,K)=D(2,I,J-1,K)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I,J-1,K)=D(3,I,J-1,K)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I,J-1,K)=D(4,I,J-1,K)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I,J-1,K)=D(5,I,J-1,K)+dreR * sav


	end do
	end do
	end do



!     k  direction      
      do k=1,NK
      do j=1,NJ1
      do i=1,NI1

         ULo=VL(3,1,i,j,k)
	   URo=VR(3,1,i,j,k)
         !!!!!!!!!!!!!!!!!
         VLo=VL(3,2,i,j,k)
	   VRo=VR(3,2,i,j,k)
         !!!!!!!!!!!!!!!!!
         WLo=VL(3,3,i,j,k)
	   WRo=VR(3,3,i,j,k)

         PLo=VL(3,4,i,j,k)
	   PRo=VR(3,4,i,j,k)

	   TLo=VL(3,5,i,j,k)
	   TRo=VR(3,5,i,j,k)
         !!!!!!!!!!!!!!!!!
        rLo=rccl(3,i,j,k)*omega(1)
        rRo=rccr(3,i,j,k)*omega(1)


          sav1=SD(3,1,I,J,K)  !aix(i,j,k)
          sav2=SD(3,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(3,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(3,I,J,K)   !aim(i,j,k)+small
		 RmiuT = VL(3,8,I,J,K)
        AlagmR = ALAGM(3,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I,J,K-1))  !by ydd


        CALL CONVECT_Roe_Common(pGam, PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR,rLo,rRo)  !output


 
	    D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I,J,K-1)=D(1,I,J,K-1)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I,J,K-1)=D(2,I,J,K-1)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I,J,K-1)=D(3,I,J,K-1)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I,J,K-1)=D(4,I,J,K-1)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I,J,K-1)=D(5,I,J,K-1)+dreR * sav

	end do
	end do
	end do

!	IF_PRECONDITION=III

	  
END SUBROUTINE CONVECT_Roe_Original





SUBROUTINE CONVECT_Roe_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR,rLo,rRo)  !output
    Use Global
    IMPLICIT NONE


REAL:: PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR  !in
REAL:: pGam                     !in  gam

REAL:: droR,drxR,dryR,drzR,dreR !out
real::rLo,rRo,rccm


REAL:: rcpcv
REAL:: DENLo,DENRo
REAL:: ga,ga1
REAL:: VAMLo,VAMRo,HLo,HRo,acl,acr
REAL:: rrorl,rrorlp1
REAL:: rm,um,vm,wm,hm,vm2,am2,am
       
REAL:: sav1n,sav2n,sav3n
REAL:: unormaln,anormaln,dacou,cef
REAL:: TTL,ubetac2
REAL:: alf,unrev,amrev,eig1,eig2,eig3,dupda,deigen
          
REAL:: dunormaln,astar,Mstar
REAL:: dW1,dU1,dU2,dU3,dP1,dW2,dW3,dW4,dW5
REAL:: dUroe,dProe

REAL:: droR1,drxR1,dryR1,drzR1,dreR1
REAL:: ulnormaln,urnormaln,rulnormaln,rurnormaln

	     rcpcv=PPF ! rcpcv=Cp-Cv=R
         DENLo=PLo/(rcpcv*TLo)
         DENRo=PRo/(rcpcv*TRo)

	   ga=pGam
	   ga1=pGam-1.0

         VAMLo=ULo*ULo+VLo*VLo+WLo*WLo   !vaml=ul*ul+vl*vl+wl*wl
         VAMRo=URo*URo+VRo*VRo+WRo*WRo   !vamr=ur*ur+vr*vr+wr*wr
         HLo=PLo*ga/(DENLo*ga1)+0.5*VAMLo-0.5*rLo**2.0 !hl=pl*ga/(rl*ga1)+0.5*vaml
         HRo=PRo*ga/(DENRo*ga1)+0.5*VAMRo-0.5*rRo**2.0 !hr=pr*ga/(rr*ga1)+0.5*vamr
         acl=sqrt(ga*PLo/DENLo)
         acr=sqrt(ga*PRo/DENRo)

        ! Roe average
          rrorl=sqrt(DENRo/DENLo)    !rrorl=sqrt(rr/rl)
          rrorlp1=1.0+rrorl
      
          rm=sqrt(DENRo*DENLo)       !rm=sqrt(rr*rl)
          um=(ULo+URo*rrorl)/rrorlp1 !um=(ul+ur*rrorl)/rrorlp1
          vm=(VLo+VRo*rrorl)/rrorlp1 !vm=(vl+vr*rrorl)/rrorlp1
          wm=(WLo+WRo*rrorl)/rrorlp1 !wm=(wl+wr*rrorl)/rrorlp1
          hm=(HLo+HRo*rrorl)/rrorlp1   !hm=(hl+hr*rrorl)/rrorlp1
          vm2=um*um+vm*vm+wm*wm      !vm2=um*um+vm*vm+wm*wm
          rccm=(rLo+rrorl*rRo)/rrorlp1          

          am2=ga1*abs(hm-0.5*vm2+0.5*rccm**2.0)    !am2=ga1*abs(hm-0.5*vm2)
          am=sqrt(am2)               !am=sqrt(am2)
!          am2i=1.0/am2                !am2i=1./am2
!          ami=1.0/am

!         surface area vectors
!

          
          sav1n=sav1/sav
          sav2n=sav2/sav
          sav3n=sav3/sav

!         engenvalues
          unormaln=(sav1n*um+sav2n*vm+sav3n*wm)  !unormal * dir not sav but savn
          anormaln=am  !*sav
          dacou=acr-acl
          cef=1.0

          IF(IF_PRECONDITION) THEN !low Ma

		 TTL= am2/ga/PPF !Tm
	     ubetac2 = 0.0
		 CALL getPrecondPara(vm2,am2,rm,TTL,RmiuT,				AlagmR,ubetac2)

			alf = (1.0-ubetac2) / 2.0
			unrev = unormaln*(1.0-alf)
			amrev = sqrt((alf*unormaln)**2.0+ubetac2*am2)

			eig1 = abs(unormaln)
			eig2 = abs(unrev+amrev)
			eig3 = abs(unrev-amrev)
			
		ELSE !high Ma limit?
			dupda=abs(sav1n*(URo-ULo)+sav2n*(VRo-VLo)+sav3n*(WRo-WLo))+abs(dacou)
		    deigen=dupda*0.5*cef  !*sav 
!			deigen=max(deigen,(abs(unormal)+abs(anormal))*0.000015)
!			eig1=abs(unormaln)
!			eig2=max(abs(unormaln+anormaln),deigen)
!			eig3=max(abs(unormaln-anormaln),deigen)

			eig1=abs(unormaln)
			eig2=abs(unormaln+anormaln)
			eig3=abs(unormaln-anormaln)
            if (eig1 <0.1) eig1 =(eig1*eig1 +0.1*0.1)/0.2
            if (eig2 <0.1) eig2 =(eig2*eig2 +0.1*0.1)/0.2
            if (eig3 <0.1) eig3 =(eig3*eig3 +0.1*0.1)/0.2
            !modified by luoj
            !if (eig1 <0.25) eig1 =(eig1*eig1 +0.25*0.25)/0.5
            !if (eig2 <0.25) eig2 =(eig2*eig2 +0.25*0.25)/0.5
            !if (eig3 <0.25) eig3 =(eig3*eig3 +0.25*0.25)/0.5


            alf=0.0
			ubetac2 = 1.0
			unrev = unormaln
			amrev = am
		ENDIF


           dunormaln=(sav1n*(URo-ULo)+sav2n*(VRo-VLo)+sav3n*(WRo-WLo)) ! dunormal=(sav1n*(ur-ul)+sav2n*(vr-vl)+sav3n*(wr-wl))
		 astar = (eig2+eig3)/2.0
		 Mstar = (eig2-eig3)/2.0/amrev

		 dW1 = DENRo-DENLo
		 dU1 = URo-ULo
		 dU2 = VRo-VLo
		 dU3 = WRo-WLo
		 dP1 = PRo-PLo
		 dW2 = rm*dU1 + dW1*um !!DENRo*URo - DENLo*ULo
		 dW3 = rm*dU2 + dW1*vm !DENRo*VRo - DENLo*VLo
		 dW4 = rm*dU3 + dW1*wm !DENRo*WRo - DENLo*WLo
		 dW5 = dP1/ga1+0.5*vm2*dW1+rm*um*dU1+rm*vm*dU2+rm*wm*dU3  
!			(PRo/ga1+0.5*DENRo*VAMRo)-
!     >			(PLo/ga1+0.5*DENLo*VAMLo)  !drE

		
		 dUroe = Mstar*dunormaln +      &
     &		(astar - (1-2.0*alf)*eig1 - alf*unormaln*Mstar)*        &
     &		(PRo-PLo)/rm/(ubetac2*am2)  !ur2 = ubetc2*am2, not amrev
		 dProe = Mstar*(PRo-PLo) +      &
     &		(astar -               eig1 + alf*unormaln*Mstar)*      &
     &		rm*dunormaln  

		 droR1 = eig1*dW1 + dUroe*rm
		 drxR1 = eig1*dW2 + dUroe*rm*um + dProe*sav1n
		 dryR1 = eig1*dW3 + dUroe*rm*vm + dProe*sav2n
		 drzR1 = eig1*dW4 + dUroe*rm*wm + dProe*sav3n
		 dreR1 = eig1*dW5 + dUroe*rm*hm + dProe*unormaln


           ulnormaln=(sav1n*ULo+sav2n*VLo+sav3n*WLo)                      ! ulnormal=(sav1*ul+sav2*vl+sav3*wl+vib(i,j,k))
           urnormaln=(sav1n*URo+sav2n*VRo+sav3n*WRo)                      ! urnormal=(sav1*ur+sav2*vr+sav3*wr+vib(i,j,k))
           rulnormaln=DENLo*ulnormaln                                   ! rulnormal=rl*ulnormal
           rurnormaln=DENRo*urnormaln                                   ! rurnormal=rr*urnormal
!           unormaln=(unormal)/sav                                     ! unormaln=(unormal-vib(i,j,k))/sav

		droR = -0.5*(rulnormaln+rurnormaln - droR1)
          drxR = -0.5*(rulnormaln*ULo+rurnormaln*URo    &
     &               +sav1n*(PLo+PRo)          - drxR1)
          dryR = -0.5*(rulnormaln*VLo+rurnormaln*VRo    &
     &               +sav2n*(PLo+PRo)          - dryR1)
          drzR = -0.5*(rulnormaln*WLo+rurnormaln*WRo    &
     &               +sav3n*(PLo+PRo)          - drzR1)
		dreR = -0.5*(rulnormaln*HLo+rurnormaln*HRo - dreR1)


END SUBROUTINE CONVECT_Roe_Common
