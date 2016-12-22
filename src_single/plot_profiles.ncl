load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

if (casein.eq."CPBL2".or.casein.eq."CPBL4") then
  plt_times=3600.*(/3.,5./)
  tavg=0.
else
  plt_times=3600.*(/2.5/)
  tavg=1800.
end if

colors=(/"blue","red","darkgreen","orange"/)
dashes=(/0,1,2,3,4/)
thickness=(/3.,1.5,1.5,1.5,1.5/)
plot=new(4,graphic)

wks   = gsn_open_wks ("eps","plot_profile")                   ; open workstation

f=addfile(filein,"r")

z=f->z
zi=f->z2
time=f->time
qt1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thl1D=new((/dimsizes(plt_times),dimsizes(z)/),float)

t_ind=new(dimsizes(plt_times),integer)
do i=0,dimsizes(t_ind)-1
  index_vec = ind(abs(plt_times(i)-time).eq.min(abs(plt_times(i)-time)))
  t_ind(i) = index_vec(0)
  delete(index_vec)
  is=ind(abs(plt_times(i)-time+tavg).eq.min(abs(plt_times(i)-time+tavg)))
  ie=ind(abs(plt_times(i)-time-tavg).eq.min(abs(plt_times(i)-time-tavg)))
  tmp=dim_avg_n(f->qt(is:ie,:),0)
  qt1D(i,:)  = dim_avg_n(f->qt(is:ie,:),0)
  thl1D(i,:) = dim_avg_n(f->thl(is:ie,:),0)
end do



; get benchmark LES data

f=addfile("/global/home/users/langhans/singlecolumn/src_ensemble/benchmark.nc","r")

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
nzles=dimsizes(zles)
thl=new((/dimsizes(Bt_ind),nzles/),float)
qt=new((/dimsizes(Bt_ind),nzles/),float)

thl=f->THL(Bt_ind,:)
qt =f->QT (Bt_ind,:)


rmse_thl=new(dimsizes(Bt_ind),float)
rmse_qt =new(dimsizes(Bt_ind),float)
rmse_thl=0.
rmse_qt=0.
do n=0,dimsizes(Bt_ind)-1
do k=0,nzles-1
  rmse_thl(n)=rmse_thl(n)+(thl(n,k)-thl1D(n,k+1))^2.
  rmse_qt(n)=rmse_qt(n)+(qt(n,k)-qt1D(n,k+1))^2.
;  print(zles(k)+" "+thl(n,k) + " "+thl1D(n,k+1)+" "+(thl(n,k)-thl1D(n,k+1))^2.)
end do
rmse_thl(n)=sqrt(rmse_thl(n)/nzles)
rmse_qt(n)=sqrt(rmse_qt(n)/nzles)
end do

if (dimsizes(Bt_ind).gt.1) then
  print("RMSE thl= "+ avg(rmse_thl))
  print("RMSE qt = "+ avg(rmse_qt))
else
  print("RMSE thl= "+ rmse_thl(0))
  print("RMSE qt = "+ rmse_qt(0))
end if

 res                   = True                       ; plot mods desired
 res@gsnFrame = False
 res@gsnDraw  = False
 res@vpHeightF = 0.5
; res@tiMainString      = "Potential temperature"             ; add title
 res@xyLineColor = "black"
 res@xyLineThicknessF = 3.
 res@xyDashPatterns = (/0,1/)
 res@vpWidthF=0.4
 res@trXMinF=299.5
 res@trXMaxF=303.
 res@trYMinF=00.
 res@trYMaxF=3000.
 res@tiXAxisString = "!F8!q!F! (K)"

 res@pmLegendDisplayMode    = "Never"            ; turn on legend
 res@pmLegendSide           = "Bottom"               ; Change location of
 res@pmLegendParallelPosF   = .7             ; move units right
 res@pmLegendWidthF         = 0.08                ; Change width and
 res@pmLegendHeightF        = 0.18                ; height of legend.
 res@lgPerimOn              = False               ; turn off box around
 res@lgLabelFontHeightF     = .015                 ; label font height

 plot(0)  = gsn_csm_xy (wks,thl,zles,res) ; create plot
 res@xyLineColor = "blue"
 plot(1)  = gsn_csm_xy (wks,thl1D,z,res) ; create plot
 overlay(plot(0),plot(1))

 res@trXMinF=0.
 res@trXMaxF=6.
 res@tiXAxisString = "q!B!t!N! (g/kg)"
 qt=qt*1000.
 qt1D=qt1D*1000.
 res@xyLineColor = "black"
 plot(2)  = gsn_csm_xy (wks,qt,zles,res) ; create plot
 res@xyLineColor = "blue"
 plot(3)  = gsn_csm_xy (wks,qt1D,z,res) ; create plot
 overlay(plot(2),plot(3))

 resP                    = True                 ; modify the panel plot
  resP@gsnFrame           = False                ; don't advance panel plot
  resP@gsnPanelLabelBar   = False                 ; add common colorbar
 resP@gsnPanelMainString = casein     ; new way of setting main title (as of NCL V6.4.0)
  gsn_panel(wks,plot((/0,2/)),(/1,2/),resP)           




