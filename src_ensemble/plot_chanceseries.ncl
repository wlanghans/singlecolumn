pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)


colors=(/"black","green","yellow","orange","red","purple","black","green","yellow","orange","red","purple"/)

improve=(/0.,0.1,0.2,0.3,0.35,0.4/)

untune=0.640

best01=0.463829

ofs=asciiread("./opti_new/of_srbf.txt",(/300,20/),"float")
ofd=asciiread("./opti_new/of_dy.txt",(/300,20/),"float")


data=new((/2*dimsizes(improve),300/),float)

do i=0,299  
ofs(i,:) = dim_min_n(ofs(0:i,:),0)
ofd(i,:) = dim_min_n(ofd(0:i,:),0)
do j=0,dimsizes(improve)-1
  data(j,i) = num(ofs(i,:).lt.(1.-improve(j))*untune)
  data(dimsizes(improve)+j,i) = num(ofd(i,:).lt.(1.-improve(j))*untune)
end do

print(i+" "+data(3,i)/20.)
end do
data=data/20.


  wks = gsn_open_wks("eps","prob_improve")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=3.
  res@xyMarkLineMode="Lines"
  res@xyLineColors= colors
  res@xyMonoDashPattern = False
  res@xyDashPatterns=(/0,0,0,0,0,0,2,2,2,2,2,2/)
  res@tiXAxisFontHeightF=0.03
  res@trYMinF = 0.

  res@tiYAxisString="Probability to improve default by X%"
  res@tiXAxisString="Number of Evaluations"


  plot  = gsn_csm_xy (wks,fspan(1.,300.,300),data(:,:),res)


  draw(plot)
  frame(wks)





