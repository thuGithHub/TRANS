SUBROUTINE COMMON_Define_BCinput_Parameter
    USE Global
    IMPLICIT NONE

    !	i_suf = 1 - 6: imin,imax,jmin,jmax,kmin,kmax
    !   i_Kind with i_Kindsub: Kind of BC
    !   i_Blk_Ex, i_Blkbd_Ex: p2p or quasip2p, correspondant domain BC
    !	  i_lcross 
    !	  i_lrev_a 
    !	  i_lrev_b 
                         
    I_Suf=1             !IJKBC(:,I_Suf),即数组第二列第一行代表边界6种位置类型：imin,imax,jmin,jmax,kmin,kmax        

    I_IBgn=2            
    I_Iend=3            
    I_JBgn=4           
    I_Jend=5           
    I_KBgn=6            
    I_Kend=7            
    
    I_Kind=8            !IJKBC(:,I_Kind),即数组第2列第8行定义为边界条件种类
    I_Kindsub=9         !IJKBC(:,I_Kindsub),即数组第2列第9行定义定义边界条件子类
    
    I_Blk_Ex=10         !IJKBC(:,I_Blk_Ex),即数组第2列第10行定义定义点对点搭接对方块号      
    I_Blkbd_Ex=11       !IJKBC(:,I_Blkbd_Ex),即数组第2列第11行定义定义点对点搭接对方块边界号
    I_Lcross=12         !IJKBC(:,I_Lcross),即数组第2列第12行定义定义本块网格是否与对应对方块是否交叉（即本块网格边界面上I方向对应对方网格块边界面上的I方向还是J方向）
    I_Lrev_A=13         !IJKBC(:,I_Lrev_A),即数组第2列第13行定义定义本网格块是否与对应对方块在I方向上相逆（即本网格块边界面上I最小值对应对方网格块I最大值还是最小值）
    I_Lrev_B=14         !IJKBC(:,I_Lrev_B),即数组第2列第14行定义定义本网格块是否与对应对方块在J方向上相逆（即本网格块边界面上I最小值对应对方网格块J最大值还是最小值）
   I_Period=15

END SUBROUTINE COMMON_Define_BCinput_Parameter
