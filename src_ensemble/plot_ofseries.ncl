pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)


colors=(/"red","black","darkgreen","red"/)

best01=0.463829
mmmodel=0.529

is=0
ie=71
ofs=asciiread("./opti_new/of_srbf.txt",(/300,20/),"float")
ofd=asciiread("./opti_new/of_dy.txt",(/300,20/),"float")


data=new((/4,300/),float)

sbtmm = False
dbtmm = False

do i=0,299  
ofs(i,:) = dim_min_n(ofs(0:i,:),0)
ofd(i,:) = dim_min_n(ofd(0:i,:),0)
data(0,i) = avg(ofs(i,:))
data(1,i) = avg(ofd(i,:))
if (.not.sbtmm.and.data(0,i).le.mmmodel) then
  print("StochRBF reaches mm skill at #eval="+i)
  sbtmm=True
end if
if (.not.dbtmm.and.data(1,i).le.mmmodel) then
  print("Dycors reaches mm skill at #eval="+i)
  dbtmm=True
end if
data(2,i) = stddev(ofs(i,:))
data(3,i) = stddev(ofd(i,:))
end do


print(min(ofs(299,:))+" "+ max(ofs(299,:)))
print(min(ofd(299,:))+" "+ max(ofd(299,:)))


  wks = gsn_open_wks("eps","trials_avgstd_dycors")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=3.
  res@xyMarkLineMode="Lines"
  res@xyLineColors= colors(1)
  res@xyMonoDashPattern = False
  res@xyDashPatterns=(/0,0,2,2/)
  res@xyMarkerColors=colors(1)
  res@xyMarker = 1
  res@xyMarkerSizeF=0.02
  res@tiXAxisFontHeightF=0.03
  res@trYMinF = 0.
  res@trYMaxF = 1.
;  res@tmYLMajorLineColor = colors(0)
;  res@tmYLMinorLineColor=colors(0)
;  res@tmYLLabelFontColor=colors(0)

;  res@tmYRLabelsOn=True
  

  res@tiYAxisString="Objective Function Value"
  res@tiXAxisString="Number of Evaluations"



  resr=res
  delete(resr@tiYAxisString)
  resr@trYMinF = -0.2
  resr@trYMaxF = 0.8
  resr@xyLineColors= colors(1)
  resr@tmYRMajorLineColor = colors(1)
  resr@tmYRMinorLineColor=colors(1)
  resr@tmYRLabelFontColor=colors(1)
  ;plot  = gsn_csm_xy2 (wks,fspan(1.,300.,300),data(0,:),data(1,:),res,resr)
  plot  = gsn_csm_xy (wks,fspan(1.,300.,300),data(1,:),res)


  ;************************************************
; add error bars
;************************************************  
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 1                          ; index for circle
  polyres@gsMarkerSizeF     = .04                        ; size 
  
; gsn_add* templates are functions that we set to dummy values. Since
; we are going to draw numerous error bars, we create two arrays to
; hold the dummy values.
  error_bar1 = new(dimsizes(data(0,:)),graphic)
  error_bar2 = new(dimsizes(data(0,:)),graphic)

 polyres@gsMarkerColor= colors(2)
 center   = gsn_add_polymarker(wks,plot,280.,0.640,polyres)
 txres               = True                            
 txres@txFontHeightF = 0.02
 txres@txFontColor= colors(2)
 text = gsn_add_text(wks,plot,"DEF",265.,0.640,txres)

; loop through the points  
  do t=0,299
; add a polymarker as the center of the error bar located at each data
; point
; the error bar is draw in two segements, above and below the true
; value of u(t). I made up these value. You would change this to 
; true error.      
   if (mod(t,2).eq.0) then
    polyres@gsMarkerColor= colors(1)
    polyres@gsLineColor = colors(1)
    ;error_bar1(t) = gsn_add_polyline(wks,plot,(/t+1.,t+1./),(/data(0,t)-0.5*data(2,t),data(0,t)+0.5*data(2,t)/),polyres)
    ;polyres@gsLineColor = colors(1)
    error_bar2(t) = gsn_add_polyline(wks,plot,(/t+1.,t+1./),(/data(1,t)-0.5*data(3,t),data(1,t)+0.5*data(3,t)/),polyres)
   end if
  end do




  draw(plot)
  frame(wks)





