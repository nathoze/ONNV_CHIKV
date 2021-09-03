
df  =data.frame(chiktiter=  data_csv$CHIKV_obs+rnorm(N)/10,
                onnvtiter = data_csv$ONNV_obs+rnorm(N)/10,
                place = factor(data_csv$location==3))

g = ggplot(data=df)+
  geom_point(aes(x=chiktiter,y=onnvtiter, color=place), size=1.1,alpha=0.4) +
  xlab('CHIKV VNT') +  ylab('ONNV VNT')+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))+ xlim(-0.5,4.5)+ylim(-0.5,4.5)+
  scale_x_continuous(breaks=seq(0,m),
                     labels=c('20', '40', '80','160','320'))+
  scale_y_continuous(breaks=seq(0,m),
                     labels=c('20', '40', '80','160','320'))+
  scale_color_manual(values=c("gray40", "brown3","brown3","gray40"))+ theme(legend.position = "none")+
  coord_equal()

print(g)