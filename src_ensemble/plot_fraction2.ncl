; ***********************************************
; xy_1.ncl
;
; Concepts illustrated:
;   - Drawing a black-and-white XY plot
;
;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************
neval=300
data1 = asciiread("p3100_dy.txt",(/neval,20/),"float")
data2 = asciiread("p3100_srbf.txt",(/neval,20/),"float")
data=new((/2,neval/10/),float)
is=0
ie=is+9
do i=0,neval/10 - 1 
data(0,i)=sum(data1(is:ie,:))
data(1,i)=sum(data2(is:ie,:))
is=is+10
ie=ie+10
end do
data=data/200.
x_vec=fspan(5,295,neval/10)
;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("eps","fracp3100")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@xyLineThicknessF = 4.0
 res@xyLineColors = (/"black","red"/)
 res@xyDashPatterns = 0
 res@trXMaxF = 300
 res@tiYAxisString    = "Fraction of trials with !F!l!BF8!t!NF! < 100 m"
 res@tiXAxisString    = "Number of evaluations"

 res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Bottom"               ; Change location of
  res@pmLegendParallelPosF   = .8                  ; move units right
  res@pmLegendOrthogonalPosF = -1.2                ; move units down
  res@pmLegendWidthF         = 0.15                ; Change width and
  res@pmLegendHeightF        = 0.18                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .03                 ; label font height
  res@xyExplicitLegendLabels = (/"Dycor","Stoch"/)         ; create explicit labels


 plot  = gsn_csm_xy (wks,x_vec,data,res) ; create plot
end

