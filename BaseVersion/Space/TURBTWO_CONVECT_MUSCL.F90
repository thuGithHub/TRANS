SUBROUTINE TURBTWO_CONVECT_MUSCL

    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    
    REAL:: KK,TKK
    REAL:: DR,DL

    INTERFACE 
        REAL FUNCTION alimiter(var1,var2)
            REAL:: var1,var2
        END FUNCTION
    END INTERFACE

    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1


!   VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm
        DO i=0,NI+1
	   DO j=0,NJ+1
	   DO k=0,NK+1
	   DO L=6,7
		  VR(1,L,i,j,k)=VL(1,L,i,j,k)
	      VR(2,L,i,j,k)=VL(2,L,i,j,k)
	      VR(3,L,i,j,k)=VL(3,L,i,j,k)
	   end do
	   end do
	   end do
	   end do

!        i direction
         Do k=1,NK1
         Do j=1,NJ1
         Do i=0,NI

!        limite TKK

         TKK=V(6,i,j,k)/V(1,i,j,k)


         DR=VL(1,6,i+1,j,k)-TKK                 ! dr=(teil(i+1,j,k)-  te(i,j,k))
         DL=TKK-VR(1,6,i,j,k)                   ! dl=(te(i,j,k)    -teir(i,j,k))

         VL(1,6,i+1,j,k)=TKK+alimiter(DR,DL)    ! teir(i,j,k)  =te(i,j,k)-alimiter(dl,dr)

	   VR(1,6,i,j,k)=TKK-alimiter(DL,DR)      ! teil(i+1,j,k)=te(i,j,k)+alimiter(dr,dl)

!        limite TOO

         TKK=V(7,i,j,k)/V(1,i,j,k)


         DR=VL(1,7,i+1,j,k)-TKK                 ! dr=(teil(i+1,j,k)-  te(i,j,k))
         DL=TKK-VR(1,7,i,j,k)                   ! dl=(te(i,j,k)    -teir(i,j,k))

         VL(1,7,i+1,j,k)=TKK+alimiter(DR,DL)    ! teir(i,j,k)  =te(i,j,k)-alimiter(dl,dr)

	   VR(1,7,i,j,k)=TKK-alimiter(DL,DR)      ! teil(i+1,j,k)=te(i,j,k)+alimiter(dr,dl)
!      

	   end do
	   end do
	   end do




!        j direction
         Do k=1,NK1
         Do j=0,NJ
         Do i=1,NI1

!        limite TKK

         TKK=V(6,i,j,k)/V(1,i,j,k)

         DR=VL(2,6,i,j+1,k)-TKK                 ! dr=(tejl(i,j+1,k)-  te(i,j,k))
         DL=TKK-VR(2,6,i,j,k)                   ! dl=(te(i,j,k)    -tejr(i,j,k))
      
         VL(2,6,i,j+1,k)=TKK+alimiter(DR,DL)    ! tejl(i,j+1,k)=te(i,j,k)+alimiter(dr,dl)
         VR(2,6,i,j,k)=TKK-alimiter(DL,DR)      ! tejr(i,j,k)  =te(i,j,k)-alimiter(dl,dr)

 
!        limite TOO

         TKK=V(7,i,j,k)/V(1,i,j,k)


         DR=VL(2,7,i,j+1,k)-TKK                 ! dr=(tejl(i,j+1,k)-  te(i,j,k))
         DL=TKK-VR(2,7,i,j,k)                   ! dl=(te(i,j,k)    -tejr(i,j,k))
      
         VL(2,7,i,j+1,k)=TKK+alimiter(DR,DL)    ! tejl(i,j+1,k)=te(i,j,k)+alimiter(dr,dl)
         VR(2,7,i,j,k)=TKK-alimiter(DL,DR)      ! tejr(i,j,k)  =te(i,j,k)-alimiter(dl,dr)

!      

	   end do
	   end do
	   end do




!        k direction
         Do k=0,NK
         Do j=1,NJ1
         Do i=1,NI1

!        limite TKK

         TKK=V(6,i,j,k)/V(1,i,j,k)

         DR=VL(3,6,i,j,k+1)-TKK                 ! dr=(tekl(i,j,k+1)-  te(i,j,k))
         DL=TKK-VR(3,6,i,j,k)                   ! dl=(te(i,j,k)    -tekr(i,j,k))
      
         VL(3,6,i,j,k+1)=TKK+alimiter(DR,DL)    ! tekl(i,j,k+1)=te(i,j,k)+alimiter(dr,dl)
         VR(3,6,i,j,k)=TKK-alimiter(DL,DR)      ! tekr(i,j,k)  =te(i,j,k)-alimiter(dl,dr)

!        limite TOO

         TKK=V(7,i,j,k)/V(1,i,j,k)

         DR=VL(3,7,i,j,k+1)-TKK                 ! dr=(tekl(i,j,k+1)-  te(i,j,k))
         DL=TKK-VR(3,7,i,j,k)                   ! dl=(te(i,j,k)    -tekr(i,j,k))
      
         VL(3,7,i,j,k+1)=TKK+alimiter(DR,DL)    ! tekl(i,j,k+1)=te(i,j,k)+alimiter(dr,dl)
         VR(3,7,i,j,k)=TKK-alimiter(DL,DR)      ! tekr(i,j,k)  =te(i,j,k)-alimiter(dl,dr)
!      

	   end do
	   end do
	   end do



END SUBROUTINE TURBTWO_CONVECT_MUSCL


