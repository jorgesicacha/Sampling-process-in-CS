library(ggplot2)
library(dplyr)
library(cowplot)
library(magick)

data <- data.frame(group=c(rep("species 1", 6), rep("species 2", 6)),
                   true_intensity = c(0,2,4,6,8,8,
                                      0,4,2,3,5,5),
                   place = c("", "mountain", "park", "forest", "road side",""))

#ggplot(data, aes(x=as.factor(place), y=true_intensity, group=as.factor(group), color=group, linetype=as.factor(group)))+
#  geom_step(direction="vh")


df = data.frame(x=rep(1:4,2),
                y=c(2,8,10,10,
                    10,4,5,7),
                group=c(rep("species 1", 4), rep("species 2", 4)),
                y1 = rep(c(0.05, 0.7, 0.2, 0.8),2),
                y2 = c(0.2, 0.8, 0.7, 0.9,
                       0.4, 0.6, 0.7, 0.8)
              )
#ggplot(df,aes(x=x,y=y))+geom_step(size=1)+scale_x_continuous(breaks=seq(0,8,2))

g1 <- ggplot(rbind(df, 
             data.frame(x = 5, y = tail((df[df$group=="species 1",])$y, 1), group="species 1", y1=0, y2=0),
             data.frame(x = 5, y = tail((df[df$group=="species 2",])$y, 1), group="species 2", y1=0, y2=0)),
       aes(x = x - 0.5, y = y, group=as.factor(group), color=as.factor(group), linetype=as.factor(group))) + 
  geom_step(size = 1)+
  scale_x_continuous(breaks = seq(0, 4, 1),
                     minor_breaks = seq(0, 4, 1) + 0.5) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line())+
  ylim(c(0,10))+
  xlab("")+
  ylab(expression(lambda[true](s)))+ 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        legend.title = element_blank())



g2 <- ggplot(rbind(df, 
                   data.frame(x = 5, y =0 , group="species 1", y1=tail((df[df$group=="species 1",])$y1, 1), y2=tail((df[df$group=="species 1",])$y2, 1)),
                   data.frame(x = 5, y = 0, group="species 2", y1=tail((df[df$group=="species 2",])$y1, 1), y2=tail((df[df$group=="species 2",])$y2, 1))),
             aes(x = x - 0.5, y = y1, group=as.factor(group), color=as.factor(group), linetype=as.factor(group))) + 
  geom_step(size = 1)+
  scale_x_continuous(breaks = seq(0, 4, 1),
                     minor_breaks = seq(0, 4, 1) + 0.5) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line())+
  ylim(c(0,1))+
  xlab("")+
  ylab(expression(zeta(s)))+ 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        legend.title = element_blank())



g3 <- ggplot(rbind(df, 
             data.frame(x = 5, y =0 , group="species 1", y1=0, y2=tail((df[df$group=="species 1",])$y2, 1)),
             data.frame(x = 5, y = 0, group="species 2", y1=0, y2=tail((df[df$group=="species 2",])$y2, 1))),
       aes(x = x - 0.5, y = y2, group=as.factor(group), color=as.factor(group), linetype=as.factor(group))) + 
  geom_step(size = 1)+
  scale_x_continuous(breaks = seq(0, 4, 1), 
                     minor_breaks = seq(0, 4, 1) + 0.5) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line())+
  ylim(c(0,1))+
  xlab("")+
  ylab(expression(psi(s)))+ 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        legend.title = element_blank())

#estimating observed intensity
d1 <- df[df$group=="species 1",]%>%
  dplyr::mutate(obs_int = y *y1*y2)%>%
  select(obs_int)

d2 <- df[df$group=="species 2",]%>%
  dplyr::mutate(obs_int = y *y1*y2)%>%
  select(obs_int)

obs_int <- t(cbind(d1,d2))

omega <- matrix(c(1,0,
                  0.5,0.5),
                nrow=2,ncol=2, byrow=T)

est_obs_int = omega %*% obs_int%>%
  t()%>%
  c()

df= cbind(df, est_obs_int)

g4 <- ggplot(rbind(df, 
                   data.frame(x = 5, y = tail((df[df$group=="species 1",])$y, 1), group="species 1", y1=0, y2=0,est_obs_int=tail((df[df$group=="species 1",])$est_obs_int, 1)),
                   data.frame(x = 5, y = tail((df[df$group=="species 2",])$y, 1), group="species 2", y1=0, y2=0, est_obs_int=tail((df[df$group=="species 2",])$est_obs_int, 1))),
             aes(x = x - 0.5, y = est_obs_int, group=as.factor(group), color=as.factor(group), linetype=as.factor(group))) + 
  geom_step(size = 1)+
  geom_step(data = rbind(data.frame(df[,-6], est_obs_int = df$y), 
                         data.frame(x = 5, y = tail((df[df$group=="species 1",])$y, 1), group="species 1", y1=0, y2=0, est_obs_int = tail((df[df$group=="species 1",])$y, 1)),
                         data.frame(x = 5, y = tail((df[df$group=="species 2",])$y, 1), group="species 2", y1=0, y2=0, est_obs_int = tail((df[df$group=="species 2",])$y, 1))),
            aes(x = x - 0.5, y = est_obs_int, group=as.factor(group), linetype=as.factor(group)), col = "grey", size = 1)+
  scale_x_continuous(breaks = seq(0, 4, 1), 
                     minor_breaks = seq(0, 4, 1) + 0.5, 
                     labels = c("","Mountains", "Park","Forest", "Road side")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line())+
  ylim(c(0,10))+
  xlab("")+
  ylab(expression(lambda[obs](s)))+
  guides(color=guide_legend(""))+ 
  theme(axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(), 
        legend.position = "none")

pimage <- axis_canvas(g4, axis = 'x') + 
  draw_image("mountains.jpeg", x = 0.5, y = - 1,scale = 4) +
  draw_image("park.jpeg", x = 1.5, y = - 1, scale = 4) +
  draw_image("forest.jpeg", x = 2.5, y = - 1, scale = 4)+
  draw_image("road_side.jpeg", x = 3.5, y = - 1, scale = 4) 

gg4 <- ggdraw(insert_xaxis_grob(g4, pimage, position = "bottom"))

figure <- ggpubr::ggarrange(g1,g2,g3,g4, 
                            nrow=4, 
                            ncol=1, 
                            common.legend = TRUE,
                            heights = c(1.8,1.8,1.8,2))

ggsave("generatingProcess.png", 
       plot = figure, 
       dpi = 1200, 
       height = 18, 
       width = 18, 
       units = "cm")



