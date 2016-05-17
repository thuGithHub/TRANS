SUBROUTINE SPACEINTERP_GetSurfCenterVars
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: KKtm

    REAL:: Den1,UUU1,VVV1,WWW1,PPP1,TTT1
    REAL:: Den,UUU,VVV,WWW,PPP,TTT    
    REAL:: alsinv
    REAL:: Tkk1,TOO1,RmiuT1
    REAL:: Tkk,TOO,RmiuT

    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
      do k=1,NK1
      do j=1,NJ1
!      do i=1,NI    !by ydd
      do i=0,NI+1
        
  	     Den1   = V(1,I-1,J,K)
	     UUU1   = V(2,I-1,J,K)/Den1
	     VVV1   = V(3,I-1,J,K)/Den1
	     WWW1   = V(4,I-1,J,K)/Den1
	     PPP1   = PP(I-1,J,K)
           TTT1   = T(I-1,j,K)

  	     Den   = V(1,I,J,K)
	     UUU   = V(2,I,J,K)/Den
	     VVV   = V(3,I,J,K)/Den
	     WWW   = V(4,I,J,K)/Den
	     PPP   = PP(I,J,K)
         TTT   = T(I,J,K)

            alsinv=1.0/(Alagm(1,i,j,k)+Alagm(1,i-1,j,k))                               !1./(alagmx(i,j,k)+alagmx(i-1,j,k))
            VL(1,1,i,j,k)=(Alagm(1,i,j,k)*UUU1+Alagm(1,i-1,j,k)*UUU )*alsinv
            VL(1,2,i,j,k)=(Alagm(1,i,j,k)*VVV1+Alagm(1,i-1,j,k)*VVV)*alsinv
            VL(1,3,i,j,k)=(Alagm(1,i,j,k)*WWW1+Alagm(1,i-1,j,k)*WWW )*alsinv
            VL(1,4,i,j,k)=(Alagm(1,i,j,k)*PPP1+Alagm(1,i-1,j,k)*PPP)*alsinv
            VL(1,5,i,j,k)=(Alagm(1,i,j,k)*TTT1+Alagm(1,i-1,j,k)*TTT )*alsinv
!        uil(i,j,k)=(alagmx(i,j,k)*vx(i-1,j,k)+alagmx(i-1,j,k)*vx(i,j,k))&*alsinv
!if(thisBlock%ID_Present_blk==4.and.k==1)then
!    write(*,*)"VL",KI,i,j,k,VL(1,1,i,j,k),UUU1,V(1,I-1,j,k),UUU,V(1,i,j,k)
!endif  
            rccl(1,i,j,k)=(Alagm(1,i,j,k)*rad(i-1,j,k)+Alagm(1,i-1,j,k)*rad(i,j,k))*alsinv
            rccr(1,i,j,k)=rccl(1,i,j,k)
        do L=1,5
            VR(1,L,i,j,k)=VL(1,L,i,j,k)
        enddo
      end do
      end do
      end do




      do k=1,NK1
!      do j=1,NJ
      do j=0,NJ+1
      do i=1,NI1

        
  	     Den1   = V(1,I,J-1,K)
	     UUU1   = V(2,I,J-1,K)/Den1
	     VVV1   = V(3,I,J-1,K)/Den1
	     WWW1   = V(4,I,J-1,K)/Den1
	     PPP1   = PP(I,J-1,K)
           TTT1   = T(I,J-1,K)

  	     Den   = V(1,I,J,K)
	     UUU   = V(2,I,J,K)/Den
	     VVV   = V(3,I,J,K)/Den
	     WWW   = V(4,I,J,K)/Den
	     PPP   = PP(I,J,K)
           TTT   = T(I,J,K)

!         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega gm muT

          alsinv=1.0/(Alagm(2,i,j,k)+Alagm(2,i,j-1,k))
          VL(2,1,i,j,k)=(Alagm(2,i,j,k)*UUU1+Alagm(2,i,j-1,k)*UUU )*alsinv
          VL(2,2,i,j,k)=(Alagm(2,i,j,k)*VVV1+Alagm(2,i,j-1,k)*VVV)*alsinv
	    VL(2,3,i,j,k)=(Alagm(2,i,j,k)*WWW1+Alagm(2,i,j-1,k)*WWW )*alsinv
	    VL(2,4,i,j,k)=(Alagm(2,i,j,k)*PPP1+Alagm(2,i,j-1,k)*PPP)*alsinv
	    VL(2,5,i,j,k)=(Alagm(2,i,j,k)*TTT1+Alagm(2,i,j-1,k)*TTT )*alsinv
!          ujl(i,j,k)=(alagmy(i,j,k)*vx(i,j-1,k)+alagmy(i,j-1,k)*vx(i,j,k))&*alsinv
        do L=1,5
            VR(2,L,i,j,k)=VL(2,L,i,j,k)
        enddo
        rccl(2,i,j,k)=(Alagm(2,i,j,k)*rad(i,j-1,k)+Alagm(2,i,j-1,k)*rad(i,j,k))*alsinv
        rccr(2,i,j,k)=rccl(2,i,j,k)
      end do
      end do
      end do
      
      !do k=1,NK
      do k=0,NK+1
      do j=1,NJ1
      do i=1,NI1

        
  	     Den1   = V(1,I,J,K-1)
	     UUU1   = V(2,I,J,K-1)/Den1
	     VVV1   = V(3,I,J,K-1)/Den1
	     WWW1   = V(4,I,J,K-1)/Den1
	     PPP1   = PP(I,J,K-1)
           TTT1   = T(I,J,K-1)

  	     Den   = V(1,I,J,K)
	     UUU   = V(2,I,J,K)/Den
	     VVV   = V(3,I,J,K)/Den
	     WWW   = V(4,I,J,K)/Den
	     PPP   = PP(I,J,K)
           TTT   = T(I,J,K)

!         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega gm muT

          alsinv=1.0/(Alagm(3,i,j,k)+Alagm(3,i,j,k-1))   !alsinv=1/(alagmz(i,j,k)+alagmz(i,j,k-1))
          VL(3,1,i,j,k)=(Alagm(3,i,j,k)*UUU1+Alagm(3,i,j,k-1)*UUU )*alsinv
          VL(3,2,i,j,k)=(Alagm(3,i,j,k)*VVV1+Alagm(3,i,j,k-1)*VVV)*alsinv
	    VL(3,3,i,j,k)=(Alagm(3,i,j,k)*WWW1+Alagm(3,i,j,k-1)*WWW )*alsinv
	    VL(3,4,i,j,k)=(Alagm(3,i,j,k)*PPP1+Alagm(3,i,j,k-1)*PPP)*alsinv
	    VL(3,5,i,j,k)=(Alagm(3,i,j,k)*TTT1+Alagm(3,i,j,k-1)*TTT )*alsinv
!          ujl(i,j,k)=(alagmy(i,j,k)*vx(i,j-1,k)+alagmy(i,j-1,k)*vx(i,j,k))&*alsinv
        do L=1,5
            VR(3,L,i,j,k)=VL(3,L,i,j,k)
        enddo
        rccl(3,i,j,k)=(Alagm(3,i,j,k)*rad(i,j,k-1)+Alagm(3,i,j,k-1)*rad(i,j,k))*alsinv
        rccr(3,i,j,k)=rccl(3,i,j,k)
      end do
      end do
      end do

!      VL(3,ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega  muT gm
!      V: r rU rV rW rH rK romega rgm
    if (IF_turb) then
      ! I direction

      do k=1,NK1
      do j=1,NJ1
      do i=0,NI+1
        
  	     Den1   = V(1,I-1,J,K)
	     Tkk1   = V(6,I-1,J,K)/Den1
	     TOO1   = V(7,I-1,J,K)/Den1
           RmiuT1 = Rmiu(i-1,j,k)

  	     Den   = V(1,I,J,K)
	     Tkk   = V(6,I,J,K)/Den
	     TOO   = V(7,I,J,K)/Den
	     RmiuT= Rmiu(i,j,k)

!         VL(3,9,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm 
          alsinv=1.0/(Alagm(1,i,j,k)+Alagm(1,i-1,j,k))                               !1./(alagmx(i,j,k)+alagmx(i-1,j,k))
          VL(1,6,i,j,k)=(Alagm(1,i,j,k)*Tkk1+Alagm(1,i-1,j,k)*Tkk )*alsinv
          VL(1,7,i,j,k)=(Alagm(1,i,j,k)*TOO1+Alagm(1,i-1,j,k)*TOO)*alsinv
          VL(1,8,i,j,k)=(Alagm(1,i,j,k)*RmiuT1+Alagm(1,i-1,j,k)*RmiuT)*alsinv
!        uil(i,j,k)=(alagmx(i,j,k)*vx(i-1,j,k)+alagmx(i-1,j,k)*vx(i,j,k))&*alsinv
        do L=6,8
            VR(1,L,i,j,k)=VL(1,L,i,j,k)
        enddo
      end do
      end do
      end do
!      VL(3,9,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega  muT gm
!      V: r rU rV rW rH rK romega rgm
      ! J direction
      do k=1,NK1
      do j=0,NJ+1
      do i=1,NI1
  	     Den1   = V(1,I,J-1,K)
	     Tkk1   = V(6,I,J-1,K)/Den1
	     TOO1   = V(7,I,J-1,K)/Den1
	     RmiuT1 = Rmiu(i,j-1,k)

  	     Den   = V(1,I,J,K)
	     Tkk   = V(6,I,J,K)/Den
	     TOO   = V(7,I,J,K)/Den
	     RmiuT= Rmiu(i,j,k)

!         VL(3,9,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega  muT gm
          alsinv=1.0/(Alagm(2,i,j,k)+Alagm(2,i,j-1,k))
          VL(2,6,i,j,k)=(Alagm(2,i,j,k)*TKK1+Alagm(2,i,j-1,k)*TKK )*alsinv
          VL(2,7,i,j,k)=(Alagm(2,i,j,k)*TOO1+Alagm(2,i,j-1,k)*TOO)*alsinv
          VL(2,8,i,j,k)=(Alagm(2,i,j,k)*RmiuT1+Alagm(2,i,j-1,k)*RmiuT)*alsinv
!          ujl(i,j,k)=(alagmy(i,j,k)*vx(i,j-1,k)+alagmy(i,j-1,k)*vx(i,j,k))&*alsinv
        do L=6,8
            VR(2,L,i,j,k)=VL(2,L,i,j,k)
        enddo
      end do
      end do
      end do

!      k direction
      do k=0,NK+1
      do j=1,NJ1
      do i=1,NI1

  	     Den1   = V(1,I,J,K-1)
	     Tkk1   = V(6,I,J,K-1)/Den1
	     TOO1   = V(7,I,J,K-1)/Den1
	     RmiuT1 = Rmiu(i,j,k-1)

  	     Den   = V(1,I,J,K)
	     Tkk   = V(6,I,J,K)/Den
	     TOO   = V(7,I,J,K)/Den
	     RmiuT= Rmiu(i,j,k)

!         VL(3,9,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm
          alsinv=1.0/(Alagm(3,i,j,k)+Alagm(3,i,j,k-1))   !alsinv=1/(alagmz(i,j,k)+alagmz(i,j,k-1))
          VL(3,6,i,j,k)=(Alagm(3,i,j,k)*Tkk1+Alagm(3,i,j,k-1)*Tkk )*alsinv
          VL(3,7,i,j,k)=(Alagm(3,i,j,k)*TOO1+Alagm(3,i,j,k-1)*TOO)*alsinv
          VL(3,8,i,j,k)=(Alagm(3,i,j,k)*RmiuT1+Alagm(3,i,j,k-1)*RmiuT )*alsinv
!          ujl(i,j,k)=(alagmy(i,j,k)*vx(i,j-1,k)+alagmy(i,j-1,k)*vx(i,j,k))&*alsinv
        do L=6,8
            VR(3,L,i,j,k)=VL(3,L,i,j,k)
        enddo
      end do
      end do
      end do
    endif

!   VR=VL 
      


END SUBROUTINE SPACEINTERP_GetSurfCenterVars
