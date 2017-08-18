load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

if (casein.eq."CPBL2".or.casein.eq."CPBL4") then
  plt_times=3600.*(/3.,5./)
  tavg=600.
else
  plt_times=3600.*(/2.5/)
  tavg=1800.
end if

colors=(/"blue"/)
dashes=(/0,1,2,3,4/)
thickness=(/3.,1.5,1.5,1.5,1.5/)
plot=new(20,graphic)

wks   = gsn_open_wks ("eps","plot_"+casein)                   ; open workstation

f=addfile(filein,"r")

z=f->z
zi=f->z2
nzm=dimsizes(z)
time=f->time
qt1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thl1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thr1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thvd1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thvd1D@_FillValue=-999.
thvd1D=0.
thgrad1D=new((/dimsizes(plt_times)/),float)

wqt1D=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wqt1D_mf=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D_mf=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wqt1D_ed=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D_ed=new((/dimsizes(plt_times),dimsizes(zi)/),float)
rhoi1D =new(dimsizes(zi),float)
BU=new((/3,dimsizes(plt_times),dimsizes(z)/),float)
WU=new((/3,dimsizes(plt_times),dimsizes(zi)/),float)

if (casein.eq."BOMEX") then
  qc1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
end if
  qtd1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
  qtd1D@_FillValue=-999.
  qtd1D=0.

n1=39
n2=24
n3=4
t_ind=new(dimsizes(plt_times),integer)
do i=0,dimsizes(t_ind)-1
  index_vec = ind(abs(plt_times(i)-time).eq.min(abs(plt_times(i)-time)))
  t_ind(i) = index_vec(0)
  is_vec=ind(abs(plt_times(i)-tavg-time).eq.min(abs(plt_times(i)+tavg-time)))
  ie_vec=ind(abs(plt_times(i)-time+tavg).eq.min(abs(plt_times(i)-time+tavg)))
  is=is_vec(0)
  ie=ie_vec(0)
  delete(index_vec)
  qt1D(i,:)  = dim_avg_n(f->qt(is:ie,:),0)
  thl1D(i,:) = dim_avg_n(f->thl(is:ie,:),0)
  thr1D(i,:) = dim_avg_n(f->thr(is:ie,:),0)
  rhoi1D(1:nzm-1)  = ((0.5*(dim_avg_n(f->rho(is:ie,0:nzm-2),0) + dim_avg_n(f->rho(is:ie,1:nzm-1),0))))
  rhoi1D(0) = 2.*rhoi1D(1) - rhoi1D(2)
  rhoi1D(nzm)= 2.*rhoi1D(nzm-1) - rhoi1D(nzm-2)
  wqt1D(i,:)  = dim_avg_n(f->qtflx(is:ie,:),0)
  wthl1D(i,:)  = dim_avg_n(f->tflx(is:ie,:),0)
  wqt1D_mf(i,:)  = 2.5e06* dim_avg_n(f->qtflx_mf(is:ie,:),0) * rhoi1D
  wthl1D_mf(i,:)  = 1005. * dim_avg_n(f->tflx_mf(is:ie,:),0) * rhoi1D
  wqt1D_ed(i,:)  = 2.5e06 * dim_avg_n(f->qtflx_ed(is:ie,:),0) * rhoi1D
  wthl1D_ed(i,:)  = 1005. * dim_avg_n(f->tflx_ed(is:ie,:),0) * rhoi1D
  BU(0,i,:)  = dim_avg_n(f->B(is:ie,n1,:),0)
  WU(0,i,:)  = dim_avg_n(f->upw(is:ie,n1,:),0)
  BU(1,i,:)  = dim_avg_n(f->B(is:ie,n2,:),0)
  WU(1,i,:)  = dim_avg_n(f->upw(is:ie,n2,:),0)
  BU(2,i,:)  = dim_avg_n(f->B(is:ie,n3,:),0)
  WU(2,i,:)  = dim_avg_n(f->upw(is:ie,n3,:),0)
  if (casein.eq."BOMEX") then
    qc1D(i,:)  = dim_avg_n(f->qcl(is:ie,:),0)
  end if

  upthv=f->upthv(is:ie,:,:)
  upqt=f->upqt(is:ie,:,:)
  fraca=f->upa(is:ie,:,:)
  fractot=dim_sum_n(fraca,1)

  nt=ie-is+1

    do k=0,nzm-1
    navg=0
    do indi=0,nt-1
      if (fractot(indi,k).ge.0.01) then
         fracsmall=0.
         fracsmallold=0.
         np2=39
         lfrac=True
         do while (np2.ge.0.and.lfrac)
          fracsmall=fracsmallold+fraca(indi,np2,k)
          if (fracsmall.gt.0.01) then
            wgt = (0.01-fracsmallold)/0.01
          else
            wgt = fraca(indi,np2,k)/0.01
          end if
          thvd1D(i,k) = thvd1D(i,k) + upthv(indi,np2,k) * wgt
          qtd1D(i,k) = qtd1D(i,k) + upqt(indi,np2,k) * wgt
          fracsmallold=fracsmall
          if (fracsmallold.gt.0.01) then
            lfrac=False
          end if
          np2=np2-1
         end do
         navg=navg+1
      end if
    end do
    if (navg.gt.0) then
       thvd1D(i,k) = thvd1D(i,k) / int2flt(navg)
       qtd1D(i,k) = qtd1D(i,k) / int2flt(navg)
    else
       thvd1D(i,k) = -999.
       qtd1D(i,k) = -999.
    end if
    end do
 
    do k=0,nzm-2
        if ((.not.ismissing(thvd1D(i,k))).and.(.not.ismissing(thvd1D(i,k+1)))) then
           tmp = 0.5*(thvd1D(i,k)+thvd1D(i,k+1))
           thvd1D(i,k)= tmp - thr1D(i,k)
           tmp = 0.5*(qtd1D(i,k)+qtd1D(i,k+1))
           qtd1D(i,k)= tmp - qt1D(i,k)
           delete(tmp)
        else
           thvd1D(i,k)= -999.
           qtd1D(i,k)= -999.
        end if
    end do
    thvd1D(i,nzm-1)= 0.
    qtd1D(i,nzm-1)= 0.
 
    delete(is_vec)
    delete(ie_vec)

    tmp=thl1D(i,:)
    zlow=1200.+i*200.
    zhigh=zlow+400.
    zlow=700.
    tmp=where(z.gt.zhigh.or.z.lt.zlow,-999.,tmp)
    tmp@_FillValue=-999.
    thgrad1D(i)=regline(doubletofloat(z),tmp)
    delete(tmp)

end do



; get benchmark LES data

f=addfile("/global/home/users/langhans/singlecolumn/src_ensemble/benchmark_SAM_BOMEXhighlarge_qt.nc","r")

if (casein.eq."CPBL2") then
   Bt_ind = (/1,2/)
end if
if (casein.eq."CPBL4") then
   Bt_ind = (/3,4/)
end if
if (casein.eq."BOMEX")
   Bt_ind = (/0/)
end if

zles=f->height1
zflx=f->height2
nzles=dimsizes(zles)
nzflx=dimsizes(zflx)
thl=new((/dimsizes(Bt_ind),nzles/),float)
qt=new((/dimsizes(Bt_ind),nzles/),float)
wqt=new((/dimsizes(Bt_ind),nzflx/),float)
wthl=new((/dimsizes(Bt_ind),nzflx/),float)
thvd=new((/dimsizes(Bt_ind),nzles/),float)
thgrad=new((/dimsizes(Bt_ind)/),float)


thl=f->THL(Bt_ind,:)
qt =f->QT (Bt_ind,:)
wthl=f->THLflx(Bt_ind,:)
wqt=f->QTflx(Bt_ind,:)
thvd=f->THVD(Bt_ind,:)
qtd=new((/dimsizes(Bt_ind),nzles/),float)
qtd=f->QTD(Bt_ind,:)
;wthl(0,0:nzflx-1)=1005. * wthl(0,0:nzflx-1) * rhoi1D(1:nzflx)
;wqt(0,0:nzflx-1)=2.5e06 * wqt(0,0:nzflx-1) * rhoi1D(1:nzflx)

if (casein.eq."BOMEX") then
  qc=new((/dimsizes(Bt_ind),nzles/),float)
  qc =f->QC(Bt_ind,:) 
end if

do i =0,dimsizes(Bt_ind)-1
tmp=thl(i,:)
zlow=1200.+i*200.
zhigh=zlow+400.
zlow=700.
tmp=where(zles.gt.zhigh.or.zles.lt.zlow,-999.,tmp)
tmp@_FillValue=-999.
thgrad(i)=regline(zles,tmp)
delete(tmp)
end do


rmse_thl=new(dimsizes(Bt_ind),float)
rmse_qt =new(dimsizes(Bt_ind),float)
rmse_thvd =new(dimsizes(Bt_ind),float)
if (casein.eq."BOMEX") then
  rmse_qc =new(dimsizes(Bt_ind),float)
  rmse_qc = 0.
  rmse_wqt =new(dimsizes(Bt_ind),float)
  rmse_wqt = 0.
  rmse_wthl =new(dimsizes(Bt_ind),float)
  rmse_wthl = 0.
end if
rmse_thl=0.
rmse_qt=0.
rmse_thvd=0.
do n=0,dimsizes(Bt_ind)-1
do k=0,nzles-2
  rmse_thl(n)=rmse_thl(n)+(thl(n,k)-thl1D(n,k))^2.
  rmse_qt(n)=rmse_qt(n)+(qt(n,k)-qt1D(n,k))^2.
  ;rmse_thvd(n)=rmse_thvd(n)+(qt(n,k)-qt1D(n,k))^2.
  if (casein.eq."BOMEX") then
    rmse_qc(n)=rmse_qc(n)+(qc(n,k)-qc1D(n,k))^2.
    rmse_wqt(n)=rmse_wqt(n)+(wqt(n,k)-wqt1D(n,k+1))^2.
    rmse_wthl(n)=rmse_wthl(n)+(wthl(n,k)-wthl1D(n,k+1))^2.
  end if
;  print(zles(k)+" "+thl(n,k) + " "+thl1D(n,k+1)+" "+(thl(n,k)-thl1D(n,k+1))^2.)
end do
rmse_thl(n)=sqrt(rmse_thl(n)/nzles)
rmse_qt(n)=sqrt(rmse_qt(n)/nzles)
if (casein.eq."BOMEX") then
  rmse_qc(n)=sqrt(rmse_qc(n)/int2flt(nzles-1))
  rmse_wqt(n)=sqrt(rmse_wqt(n)/int2flt(nzles-1))
  rmse_wthl(n)=sqrt(rmse_wthl(n)/int2flt(nzles-1))
end if
end do

;if (dimsizes(Bt_ind).gt.1) then
;  print("RMSE thl= "+ avg(rmse_thl))
;  print("RMSE qt = "+ avg(rmse_qt))
;else
;  print("RMSE thl= "+ rmse_thl(0))
;  print("RMSE qt = "+ rmse_qt(0))
;  print("RMSE qc = "+ rmse_qc(0))
;end if

 res                   = True                       ; plot mods desired
 res@gsnFrame = False
 res@gsnDraw  = False
 res@vpHeightF = 0.5
; res@tiMainString      = "Potential temperature"             ; add title
 res@xyLineColor = "black"
 res@xyLineThicknessF = 3.
 res@xyDashPattern = 0
 res@vpWidthF=0.6
 res@trYMinF=00.
 if (casein.eq."BOMEX") then
   res@trYMaxF=2000.
 else
   res@trYMaxF=3000.
 end if
 res@tiYAxisString = "Height (m)"
 res@tiXAxisString = "!F8!q!F! (K)"

 res@pmLegendDisplayMode    = "Never"            ; turn on legend
 res@pmLegendSide           = "Bottom"               ; Change location of
 res@pmLegendParallelPosF   = .7             ; move units right
 res@pmLegendWidthF         = 0.08                ; Change width and
 res@pmLegendHeightF        = 0.18                ; height of legend.
 res@lgPerimOn              = False               ; turn off box around
 res@lgLabelFontHeightF     = .015                 ; label font height


 res@tmXBLabelsOn         = True            ; do not draw bottom labels
 res@tiXAxisOn            = True
 if (casein.eq."BOMEX") then
 res@trXMinF=298.5
 ;res@trXMaxF=312.
 res@trXMaxF=308.
 else
 res@trXMinF=300.
 res@trXMaxF=303.
 
 grad_txt_les=sprintf("%5.3f",avg(thgrad)*1000.) +" K/km"
 grad_txt_1D=sprintf("%5.3f",avg(thgrad1D)*1000.) +" K/km"

 end if
 plot(0)  = gsn_csm_xy (wks,thl,zles,res) ; create plot
 res@xyLineColor = colors(0)

 res@tmYLLabelsOn         = False              ; do not draw bottom labels
 res@tiYAxisOn            = False

 plot(1)  = gsn_csm_xy (wks,thl1D,z,res) ; create plot


 if (casein.ne."BOMEX") then
   txres               = True                            
   txres@txFontHeightF = 0.03           ; Set the font height
   txres@txFontColor = "black"
   text = gsn_add_text(wks,plot(0),grad_txt_les,302., \
                                    1500. ,txres)
   txres@txFontColor = colors(0)
   text = gsn_add_text(wks,plot(0),grad_txt_1D,302., \
                                    1200. ,txres)
 end if

 overlay(plot(0),plot(1))

 if (casein.eq."BOMEX") then
 res@trXMinF=9.
 res@trXMaxF=17.5
 ;res@trXMinF=14.5
 else
 res@trXMinF=2.
 res@trXMaxF=6.
 end if
 res@tiXAxisString = "q!B!t!N! (g/kg)"
 qt=qt*1000.
 qt1D=qt1D*1000.
 res@xyLineColor = "black"
 plot(2)  = gsn_csm_xy (wks,qt,zles,res) ; create plot
 res@xyLineColor = colors(0)
 plot(3)  = gsn_csm_xy (wks,qt1D,z,res) ; create plot
 overlay(plot(2),plot(3))

 if (casein.eq."BOMEX") then
 res@trXMinF=0.
 res@trXMaxF=1.e-2
 res@tiXAxisString = "q!B!c!N! (g/kg)"
 qc=qc*1000.
 qc1D=qc1D*1000.
 res@xyLineColor = "black"
 plot(4)  = gsn_csm_xy (wks,qc,zles,res) ; create plot
 res@xyLineColor = colors(0)
 plot(5)  = gsn_csm_xy (wks,qc1D,z,res) ; create plot
 overlay(plot(4),plot(5))
 end if

 if (casein.eq."BOMEX") then
 res@trXMinF=-40.
 res@trXMaxF=20.
 ;res@trXMinF=14.5
 else
 res@trXMinF=-40.
 res@trXMaxF=100.
 end if
 res@tiXAxisString = "<w !F8!q!F!!B!l!N!> (W m!S!-2!N!)"
 res@xyLineColor = "black"
 res@gsnXRefLine = 0.0
 plot(6)  = gsn_csm_xy (wks,wthl,zflx,res) ; create plot
 res@xyLineColor = colors(0)
 res@xyDashPattern = 0
 plot(7)  = gsn_csm_xy (wks,wthl1D,zi,res) ; create plot
 overlay(plot(6),plot(7))
 res@xyDashPattern = 0
 if (casein.ne."BOMEX") then
 res@tmYLLabelsOn         = True              ; do not draw bottom labels
 res@tiYAxisOn            = True
 end if
 plot(14)  = gsn_csm_xy (wks,wthl1D_ed,zi,res) ; create plot
 res@xyDashPattern = 14
 plot(10)  = gsn_csm_xy (wks,wthl1D_mf,zi,res) ; create plot
 if (casein.ne."BOMEX") then
 res@tmYLLabelsOn         = False              ; do not draw bottom labels
 res@tiYAxisOn            = False
 end if
 overlay(plot(14),plot(10))
 res@xyDashPattern = 0
 if (casein.eq."BOMEX") then
 res@trXMinF=-20.
 res@trXMaxF=160.
 ;res@trXMinF=14.5
 else
 res@trXMinF=-20.
 res@trXMaxF=90.
 end if
 res@tiXAxisString = "<w q!B!t!N!> (W m!S!-2!N!)"
 res@xyLineColor = "black"
 res@xyDashPattern = 0
 if (casein.eq."BOMEX") then
 res@tmYLLabelsOn         = True              ; do not draw bottom labels
 res@tiYAxisOn            = True
 end if
 plot(8)  = gsn_csm_xy (wks,wqt,zflx,res) ; create plot
 res@xyLineColor = colors(0)
 plot(9)  = gsn_csm_xy (wks,wqt1D,zi,res) ; create plot
 overlay(plot(8),plot(9))
 res@tmYLLabelsOn         = False              ; do not draw bottom labels
 res@tiYAxisOn            =  False
 res@tiXAxisOn            =  True
 plot(15)  = gsn_csm_xy (wks,wqt1D_ed,zi,res) ; create plot
 res@xyDashPattern = 14
 plot(11)  = gsn_csm_xy (wks,wqt1D_mf,zi,res) ; create plot
 overlay(plot(15),plot(11))

 res@trXMinF=0.
 res@trXMaxF=3.
 res@xyMonoDashPattern = False
 ;res@xyDashPatterns = (/0,1,2/)
 ;res@tiXAxisString = "w!B!u!N! (m s!S!-1!N!)"
 ;plot(12)  = gsn_csm_xy (wks,WU(:,0,:),zi,res) ; create plot
 if (casein.eq."BOMEX") then
 res@trXMinF=-1.0
 res@trXMaxF=1.0
 else
 res@trXMinF=-1.0
 res@trXMaxF=1.0
 end if
 res@tiXAxisString = "!F8!q!F!!B!v!N!' (K)"
 res@xyDashPattern =0 
 res@xyMonoDashPattern = True
 res@xyLineColor = "black"
 plot(13)  = gsn_csm_xy (wks,thvd,zles,res) ; create plot
 res@xyLineColor = colors(0)
 plot(16)  = gsn_csm_xy (wks,thvd1D,z,res) ; create plot
 overlay(plot(13),plot(16))

 if (casein.eq."BOMEX") then
 res@trXMinF=0.
 res@trXMaxF=2.
 else
 res@trXMinF=0.
 res@trXMaxF=1.
 end if
 res@xyMonoDashPattern = False
 res@tiXAxisString = "q!B!t!N!' (g/kg)"
 res@xyDashPattern =0 
 res@xyMonoDashPattern = True
 res@xyLineColor = "black"
 qtd=qtd*1000.
 qtd1D=qtd1D*1000.
 plot(17)  = gsn_csm_xy (wks,qtd,zles,res) ; create plot
 res@xyLineColor = colors(0)
 plot(18)  = gsn_csm_xy (wks,qtd1D,z,res) ; create plot
 overlay(plot(17),plot(18))

 resP                    = True                 ; modify the panel plot
  resP@gsnFrame           = False                ; don't advance panel plot
  resP@gsnPanelLabelBar   = False                 ; add common colorbar
 resP@gsnPanelMainString = casein     ; new way of setting main title (as of NCL V6.4.0)
  if (casein.eq."BOMEX") then
  gsn_panel(wks,plot((/0,2,4,6,8,14,15,13,17/)),(/3,4/),resP)           
;  gsn_panel(wks,plot((/0,6,12,13/)),(/1,4/),resP)           
  else
  gsn_panel(wks,plot((/0,2,6,8,14,15,13,17/)),(/2,4/),resP)           
  end if




